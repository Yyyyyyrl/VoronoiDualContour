import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

def visualize_points_from_csv(filename):
    points = []

    # Read the points from the CSV file
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            x, y, z = float(row['x']), float(row['y']), float(row['z'])
            points.append((x, y, z))

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the points in red
    xs, ys, zs = zip(*points)
    ax.scatter(xs, ys, zs, color='red', label='Points')

    # For each point, plot a small cube in a different color
    for i, (x, y, z) in enumerate(points):
        # Define cube size
        size = 0.5
        # Define cube vertices
        r = [-size / 2, size / 2]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s - e)) == r[1] - r[0]:
                ax.plot3D(*zip(s + np.array([x, y, z]), e + np.array([x, y, z])), color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title("3D Points and Corresponding Cubes")
    plt.show()

# Call the function with the CSV file name
visualize_points_from_csv('fuel-crop.csv')