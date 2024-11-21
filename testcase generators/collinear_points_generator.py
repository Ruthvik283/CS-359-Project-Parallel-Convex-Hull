import numpy as np

# Collinear points generation
def generate_collinear_points(n, slope, intercept, x_range):
    x = np.linspace(x_range[0], x_range[1], n)  # Generate `n` evenly spaced x-coordinates within the range
    y = slope * x + intercept  # Compute corresponding y-coordinates using the line equation y = slope * x + intercept
    points = np.column_stack((x, y)).astype(int)  # Combine x and y into an array of integer points
    return points

# Save points to file in the input format for convex hull
def save_points_to_input_format(points, filename):
    with open(filename, 'w') as f:
        # f.write(f"1\n")  # 1 test case
        f.write(f"{len(points)}\n")  # number of points
        for point in points:
            f.write(f"{point[0]} {point[1]}\n")  # write each point as integers

# Usage
if __name__ == "__main__":
    # Generate 100000 collinear integer points
    n = int(1e5)
    slope = 3  # Example slope
    intercept = 5  # Example intercept
    x_range = (-10000, 1000)  # X coordinate range

    # Generate collinear points
    collinear_points = generate_collinear_points(n, slope, intercept, x_range)
    save_points_to_input_format(collinear_points, "collinear_points_input.txt")
    print("Generated input file with 100000 collinear integer points.")
