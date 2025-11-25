import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import ceil, sqrt

pairs = [
    # Pair 1
    (
        [(100.535, 96.438, 160.045),
         (99.7385, 97.5, 159.356),
         (100.234, 96.8413, 159.75)],
        [(100.5, 96.5, 160.25),
         (99.75, 97.5, 159.5),
         (100.5, 96.5, 159.5)]
    ),
    # Pair 2
    (
        [(109.808, 93.9357, 160.908),
         (109.769, 93.8929, 161.435),
         (109.761, 93.863, 161.44)],
        [(109.75, 93.75, 160.75),
         (109.5, 94.25, 161.5),
         (109.5, 93.5, 161.5)]
    ),
    # Pair 3
    (
        [(143.616, 136.457, 68.4675),
         (144.917, 135.949, 69.0509),
         (144.226, 136.185, 68.775)],
        [(143.5, 136.5, 68.5),
         (144.75, 135.75, 69.25),
         (144.5, 136.5, 68.5)]
    )
    # Add more pairs as needed
]

def plot_triangle(ax, tri):
    x, y, z = zip(*tri)
    ax.plot(list(x) + [x[0]],
            list(y) + [y[0]],
            list(z) + [z[0]],
            marker='o')

def connect_corresponding(ax, tri_a, tri_b):
    for p, q in zip(tri_a, tri_b):
        ax.plot([p[0], q[0]],
                [p[1], q[1]],
                [p[2], q[2]],
                linestyle='--')

num_pairs = len(pairs)
cols = int(ceil(sqrt(num_pairs)))
rows = int(ceil(num_pairs / cols))
fig = plt.figure(figsize=(6 * cols, 6 * rows))

for i, (prim, dual) in enumerate(pairs, start=1):
    ax = fig.add_subplot(rows, cols, i, projection='3d')
    plot_triangle(ax, prim)
    plot_triangle(ax, dual)
    connect_corresponding(ax, prim, dual)
    ax.set_title(f"Triangle Pair {i}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

plt.tight_layout()
plt.show()
