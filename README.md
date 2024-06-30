# Graph Visualizer

This project visualizes graphs by taking edge lists as input and generating a BMP image. The implementation includes functionalities for reading graph data, processing graph layouts, and generating BMP images.

## Build

   ```bash
   cd <root project directory>
   mkdir build
   cd build
   cmake ..
   make
   ```

## Usage

1. **Prepare your input file:**

   Create a `.txt` file in the `examples` folder. The input file should contain the number of vertices and edges on the first line, followed by the edges, one per line. For example:

   ```
   4 4
   1 2
   2 3
   3 4
   4 1
   ```

2. **Run the program:**

   ```sh
   ./lab4
   ```

3. **View the results**

   Find the results in the folder `output`.

### Description of Core Components

- **bmpHandler.cc:** Handles BMP file reading and writing.
- **graphHandler.cc:** Manages the graph data structures and algorithms for layout processing.
- **graphMath.cc:** Contains mathematical functions and algorithms for graph layout optimization.
- **graphVisualizator.cc:** Handles the visualization of the graph and drawing vertices and edges.

## Examples

Some examples can be found in the folder `output`.

## References

The algorithm used in this project is described in 'A Fast Multi-Scale Method for Drawing Large Graphs' by David Harel and Yehuda Koren.
