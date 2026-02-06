
import numpy as np
import matplotlib.pyplot as plt


N = 8
basis_image = np.zeros((N*N, N*N))

for u in range(N):
    for v in range(N):
        for x in range(N):
            for y in range(N):
                value = np.exp(-1j * 2 * np.pi * ((u*x + v*y) / N))
                basis_image[u*N + x, v*N + y] = np.real(value)

plt.imshow(basis_image, cmap='gray')
plt.title("8×8 2-D DFT Basis (64×64)")
plt.colorbar()
plt.show()

img = np.zeros((64, 64))

x0 = int(input("Enter top-left x: "))
y0 = int(input("Enter top-left y: "))
w  = int(input("Enter width: "))
h  = int(input("Enter height: "))

img[x0:x0+h, y0:y0+w] = 1

plt.imshow(img, cmap='gray')
plt.title("Binary Rectangle Image")
plt.show()

def dft2d(image):
    M, N = image.shape
    F = np.zeros((M, N), dtype=complex)

    for u in range(M):
        for v in range(N):
            for x in range(M):
                for y in range(N):
                    F[u, v] += image[x, y] * np.exp(
                        -1j * 2 * np.pi * ((u*x/M) + (v*y/N))
                    )
    return F

F = dft2d(img)

plt.imshow(np.log(1 + np.abs(F)), cmap='gray')
plt.title("2-D DFT Magnitude (Uncentered)")
plt.colorbar()
plt.show()


centered_img = np.zeros_like(img)

for x in range(64):
    for y in range(64):
        centered_img[x, y] = img[x, y] * ((-1) ** (x + y))

F_centered = dft2d(centered_img)

plt.imshow(np.log(1 + np.abs(F_centered)), cmap='gray')
plt.title("2-D DFT Magnitude (Centered)")
plt.colorbar()
plt.show()
