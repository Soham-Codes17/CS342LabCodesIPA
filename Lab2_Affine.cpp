
#include <bits/stdc++.h>
using namespace std;

const double PI = 3.14159265358979323846;

struct Image {
    int w, h, maxv;
    vector<unsigned char> data;
};

string readToken(istream &in) {
    string tok;
    char ch;
    while (in.get(ch)) {
        if (!isspace((unsigned char)ch)) {
            in.unget();
            break;
        }
    }
    while (in.peek() == '#') {
        string line;
        getline(in, line);
        while (in.get(ch)) {
            if (!isspace((unsigned char)ch)) {
                in.unget();
                break;
            }
        }
    }
    while (in.get(ch)) {
        if (isspace((unsigned char)ch)) break;
        tok.push_back(ch);
    }
    return tok;
}

Image readPPM(const string &filename) {
    ifstream in(filename, ios::binary);
    if (!in) throw runtime_error("Cannot open input file");

    string magic;
    in >> magic;
    if (magic != "P6") throw runtime_error("Only P6 PPM supported");

    int w = stoi(readToken(in));
    int h = stoi(readToken(in));
    int maxv = stoi(readToken(in));
    in.get();

    Image img{w, h, maxv, vector<unsigned char>(3 * w * h)};
    in.read(reinterpret_cast<char*>(img.data.data()), img.data.size());
    return img;
}

void writePPM(const string &filename, const Image &img) {
    ofstream out(filename, ios::binary);
    out << "P6\n" << img.w << " " << img.h << "\n" << img.maxv << "\n";
    out.write(reinterpret_cast<const char*>(img.data.data()), img.data.size());
}


using Mat = array<array<double,3>,3>;

Mat matMul(const Mat &A, const Mat &B) {
    Mat C{};
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

array<double,3> matVec(const Mat &A, const array<double,3> &v) {
    return {
        A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2],
        A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2],
        A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2]
    };
}

double det3(const Mat &m) {
    return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
         - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
         + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}

Mat inverseMat(const Mat &m, bool &ok) {
    Mat inv{};
    double d = det3(m);
    if (fabs(d) < 1e-12) { ok = false; return inv; }
    ok = true;

    inv[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/d;
    inv[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])/d;
    inv[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])/d;

    inv[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])/d;
    inv[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/d;
    inv[1][2]=(m[0][2]*m[1][0]-m[0][0]*m[1][2])/d;

    inv[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/d;
    inv[2][1]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])/d;
    inv[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])/d;

    return inv;
}

array<double,3> bilinearSample(const Image &img, double x, double y) {
    int x0 = floor(x), y0 = floor(y);
    int x1 = x0 + 1, y1 = y0 + 1;
    double wx = x - x0, wy = y - y0;

    auto get = [&](int X, int Y, int c) {
        if (X < 0 || Y < 0 || X >= img.w || Y >= img.h) return 0.0;
        return (double)img.data[(Y * img.w + X) * 3 + c];
    };

    array<double,3> out{};
    for (int c = 0; c < 3; c++) {
        double v00 = get(x0, y0, c);
        double v10 = get(x1, y0, c);
        double v01 = get(x0, y1, c);
        double v11 = get(x1, y1, c);
        out[c] = (v00*(1-wx) + v10*wx)*(1-wy)
               + (v01*(1-wx) + v11*wx)*wy;
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input.ppm output.ppm\n";
        return 1;
    }

    Image src = readPPM(argv[1]);

    double sx, sy, theta_deg, tx, ty, shx, shy;
    cout << "Enter sx sy rotation(deg) tx ty shx shy:\n";
    cin >> sx >> sy >> theta_deg >> tx >> ty >> shx >> shy;

    double theta = theta_deg * PI / 180.0;

    Mat S = {{{sx,0,0},{0,sy,0},{0,0,1}}};
    Mat Sh = {{{1,shx,0},{shy,1,0},{0,0,1}}};
    Mat R = {{{cos(theta),-sin(theta),0},{sin(theta),cos(theta),0},{0,0,1}}};
    Mat T = {{{1,0,tx},{0,1,ty},{0,0,1}}};

    Mat M = matMul(T, matMul(R, matMul(Sh, S)));

    vector<array<double,3>> corners = {
        {0,0,1},{(double)src.w-1,0,1},
        {0,(double)src.h-1,1},{(double)src.w-1,(double)src.h-1,1}
    };

    double minX=1e9,minY=1e9,maxX=-1e9,maxY=-1e9;
    for (auto &p : corners) {
        auto t = matVec(M, p);
        minX = min(minX, t[0]);
        minY = min(minY, t[1]);
        maxX = max(maxX, t[0]);
        maxY = max(maxY, t[1]);
    }

    int outW = ceil(maxX - minX) + 1;
    int outH = ceil(maxY - minY) + 1;

    Image out{outW, outH, src.maxv, vector<unsigned char>(3*outW*outH, 0)};

    bool ok;
    Mat Minv = inverseMat(M, ok);
    if (!ok) {
        cerr << "Transformation matrix not invertible\n";
        return 1;
    }

    for (int j = 0; j < outH; j++) {
        for (int i = 0; i < outW; i++) {
            auto p = matVec(Minv, {minX + i, minY + j, 1});
            double x = p[0], y = p[1];

            if (x >= 0 && y >= 0 && x < src.w && y < src.h) {
                auto c = bilinearSample(src, x, y);
                int idx = (j * outW + i) * 3;

                for (int k = 0; k < 3; k++) {
                    double val = c[k];
                    if (val < 0) val = 0;
                    if (val > 255) val = 255;
                    out.data[idx + k] = (unsigned char)val;
                }
            }
        }
    }

    writePPM(argv[2], out);
    cout << "Done. Output written to " << argv[2] << "\n";
    return 0;
}

