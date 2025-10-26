import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def relu(x, offset=0):
    """ReLU function: max(0, x - offset)"""
    return np.maximum(0, x - offset)

def convolve_relus(offsets, x_range=(-5, 10), num_points=1000):
    """
    Convolve N ReLU functions with specified offsets.
    
    Parameters:
    -----------
    offsets : list or array
        List of offset values for each ReLU function
    x_range : tuple
        Range of x values (min, max)
    num_points : int
        Number of points for discretization
    
    Returns:
    --------
    x_vals : array
        X values for individual ReLU functions
    relu_funcs : list of arrays
        Individual ReLU function values
    conv_x : array
        X values for convolution result
    conv_result : array
        Convolution result values
    """
    # Generate x values
    x_vals = np.linspace(x_range[0], x_range[1], num_points)
    dx = x_vals[1] - x_vals[0]
    
    # Calculate individual ReLU functions
    relu_funcs = [relu(x_vals, offset) for offset in offsets]
    
    # Perform sequential convolution
    conv_result = relu_funcs[0].copy()
    
    for i in range(1, len(relu_funcs)):
        # Discrete convolution with proper scaling
        conv_result = signal.convolve(conv_result, relu_funcs[i], mode='full') * dx
    
    # Generate x values for convolution result
    conv_length = len(conv_result)
    conv_x_range = (x_range[0] * len(offsets), x_range[1] * len(offsets))
    conv_x = np.linspace(conv_x_range[0], conv_x_range[1], conv_length)
    
    return x_vals, relu_funcs, conv_x, conv_result

def plot_relu_convolution(offsets, x_range=(-5, 10), num_points=1000):
    """
    Plot N ReLU functions and their convolution.
    
    Parameters:
    -----------
    offsets : list or array
        List of offset values for each ReLU function
    x_range : tuple
        Range of x values (min, max)
    num_points : int
        Number of points for discretization
    """
    x_vals, relu_funcs, conv_x, conv_result = convolve_relus(offsets, x_range, num_points)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot individual ReLU functions
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(offsets)))
    for i, (relu_func, offset) in enumerate(zip(relu_funcs, offsets)):
        ax1.plot(x_vals, relu_func, label=f'ReLU (offset={offset})', 
                color=colors[i], linewidth=2)
    
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('f(x)', fontsize=12)
    ax1.set_title(f'Individual ReLU Functions (N={len(offsets)})', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Plot convolution result
    ax2.plot(conv_x, conv_result, color='crimson', linewidth=2.5, label='Convolution Result')
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('Amplitude', fontsize=12)
    ax2.set_title('Convolution of All ReLU Functions', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # Print statistics
    print(f"Number of ReLU functions: {len(offsets)}")
    print(f"Offsets: {offsets}")
    print(f"Convolution maximum: {np.max(conv_result):.4f}")
    print(f"Convolution result length: {len(conv_result)} points")

# Example usage
if __name__ == "__main__":
    # Example 1: Equidistant offsets
    print("Example 1: Equidistant offsets")
    offsets_equi = [-2, -1, 0, 1, 2]
    plot_relu_convolution(offsets_equi)
    
    # Example 2: Non-equidistant offsets
    print("\nExample 2: Non-equidistant offsets")
    offsets_non_equi = [-2.5, -0.8, 0.3, 1.7, 3.2]
    plot_relu_convolution(offsets_non_equi)
    
    # Example 3: More ReLU functions
    print("\nExample 3: 8 ReLU functions")
    offsets_many = np.linspace(-3, 3, 8)
    plot_relu_convolution(offsets_many)
    
    # Example 4: Custom offsets
    print("\nExample 4: Custom random offsets")
    np.random.seed(42)
    offsets_random = sorted(np.random.uniform(-4, 4, 10))
    plot_relu_convolution(offsets_random, x_range=(-6, 12))