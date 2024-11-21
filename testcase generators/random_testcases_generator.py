import numpy as np

# Random integer points generation
def generate_random_integer_points(n, x_range, y_range):
    return np.random.randint(low=x_range[0], high=x_range[1] + 1, size=(n, 2))

# Save points to file in the input format for convex hull
def save_points_to_input_format(points, filename):
    with open(filename, 'w') as f:
        # f.write(f"1\n")  # 1 test case
        f.write(f"{len(points)}\n")  # number of points
        for point in points:
            f.write(f"{point[0]} {point[1]}\n")  # write each point as integers

# Usage
if __name__ == "__main__":
    # Generate 100000 random integer points

    n = int(5e6)
    x_range = (-1000000, 1000000)
    y_range = (-1000000, 1000000)

    s="input"+str(n)+".txt"

    random_integer_points = generate_random_integer_points(n, x_range, y_range)
    save_points_to_input_format(random_integer_points, s)
    print(f"Generated input file with {n} random integer points.")
