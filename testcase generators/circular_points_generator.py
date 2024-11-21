import numpy as np

# Circular integer points generation
def generate_circular_integer_points(n, radius, center, seed=None):
    if seed is not None:
        np.random.seed(seed)  # Set the random seed for reproducibility

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)
    points = np.column_stack((x, y))
    return np.round(points).astype(int)  # Round to integers

# Save points to file in the input format for convex hull
def save_points_to_input_format(points, filename):
    with open(filename, 'w') as f:
        # f.write(f"1\n")  # 1 test case
        f.write(f"{len(points)}\n")  # number of points
        for point in points:
            f.write(f"{point[0]} {point[1]}\n")  # write each point as integers

# Main usage
if __name__ == "__main__":
    n = int(5e5)  # Number of points
    seed = 42  # Seed for reproducibility

    # Generate circular points
    radius = 100000
    center = (0, 0)
    circular_points = generate_circular_integer_points(n, radius, center, seed)
    save_points_to_input_format(circular_points, "circular_points_input.txt")
    print("Generated circular points and saved to 'circular_points_input.txt'.")
