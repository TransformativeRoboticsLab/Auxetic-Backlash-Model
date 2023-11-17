import math


class SeriesManipulator:
    def __init__(self, num_links):
        self.num_links = num_links
        if self.num_links < 2:
            raise ValueError("The number of links should be at least 2.")
        self.angles = [0] * (num_links - 1)  # One less angle, as the last angle will be determined automatically

    def set_joint_angles(self, angles):
        """Set angles for each joint, excluding the last joint."""
        if len(angles) != self.num_links - 1:
            raise ValueError(f"Expected {self.num_links - 1} angles, got {len(angles)}")
        self.angles = angles

    def get_link_end_positions(self):
        """Get the end position of each link."""
        positions = [(0, 0)]  # Starting at the origin
        x, y = 0, 0
        for angle in self.angles:
            x += math.cos(angle)
            y += math.sin(angle)
            positions.append((x, y))

        # Calculate the angle for the last link to reach (5,0)
        dx = 5 - x
        dy = 0 - y
        last_angle = math.atan2(dy, dx)
        x += math.cos(last_angle)
        y += math.sin(last_angle)
        positions.append((x, y))

        return positions


if __name__ == "__main__":
    # Example usage
    manipulator = SeriesManipulator(num_links=10)
    # For this example, I'm setting 0 rad angles. You'll need to figure out appropriate angles to get to (5,0)
    manipulator.set_joint_angles(
        [math.radians(0) for _ in range(9)])  # Setting each angle to 0 degrees as a placeholder
    positions = manipulator.get_link_end_positions()

    for i, (x, y) in enumerate(positions):
        print(f"End position of link {i}: (x={x:.2f}, y={y:.2f})")
