import edu.princeton.cs.algs4.StdDraw;
import edu.princeton.cs.algs4.StdRandom;
import edu.princeton.cs.algs4.WeightedQuickUnionUF;

public class Percolation {
    private boolean[][] grid; // Grid to represent open/closed sites
    private int gridSize; // Size of the grid
    private WeightedQuickUnionUF wquFind; // Data structure to track connected components
    private int virtualTop; // Virtual top site for connectivity checking
    private int virtualBottom; // Virtual bottom site for connectivity checking

/**
 * Constructor for the Percolation class.
 * Initializes the percolation grid and sets up the union-find data structure to manage connectivity.
 *
 * @param n The size of the percolation grid. The grid will be n-by-n.
 * @throws IllegalArgumentException if n is less than or equal to 0.
 */
public Percolation(int n) {
    // Check if the given grid size is valid. If n is less than or equal to 0, throw an exception.
    if (n <= 0) {
        throw new IllegalArgumentException("Grid size must be greater than 0.");
    }

    // Set the grid size based on the provided parameter.
    gridSize = n;

    // Initialize the grid. Each element in the grid represents a site.
    // By default, all sites are initialized to false, indicating that they are closed.
    grid = new boolean[n][n];

    // Initialize the union-find data structure with a size to accommodate all sites in the grid
    // plus two additional sites for the virtual top and bottom.
    wquFind = new WeightedQuickUnionUF(n * n + 2); // +2 for virtual top and bottom

    // Calculate the indices for the virtual top and bottom sites in the union-find structure.
    // These virtual sites are used to simplify the check for percolation.
    virtualTop = n * n;       // Index for the virtual top site
    virtualBottom = n * n + 1; // Index for the virtual bottom site
}


/**
 * Initializes connections between virtual top and bottom sites and the top and bottom rows of the grid.
 * This method is called after constructing the Percolation object to set up the initial state of the grid.
 */
private void initializeConnections() {
    for (int i = 0; i < gridSize; i++) {
        // Connect each site in the top row to the virtual top site.
        // This is done by calling union on the index of the virtual top and each top row site.
        wquFind.union(virtualTop, i);

        // Connect each site in the bottom row to the virtual bottom site.
        // The xyTo1D method converts the 2D coordinate (bottom row, ith column) to a 1D index for the union-find structure.
        // The bottom row is represented by 'gridSize - 1' since array indices start at 0.
        wquFind.union(virtualBottom, xyTo1D(gridSize - 1, i));
    }
}

    

/**
 * Converts 2D grid coordinates to a 1D array index.
 * This is used to map the 2D grid coordinates (row and column) to the 1D structure
 * of the union-find data structure.
 *
 * @param row The row index in the grid.
 * @param col The column index in the grid.
 * @return The corresponding index in a 1D array representation of the grid.
 */
private int xyTo1D(int row, int col) {
    // The 1D index is calculated as 'row' times the size of the grid plus 'col'.
    // This ensures a unique index for each (row, col) pair in the 2D grid.
    return row * gridSize + col;
}

/**
 * Validates that the specified coordinates are within the bounds of the grid.
 * This method ensures that any access to the grid using the provided row and column
 * indices is safe and will not cause an ArrayIndexOutOfBoundsException.
 *
 * @param row The row index to validate.
 * @param col The column index to validate.
 * @throws IndexOutOfBoundsException if the row or column is outside the grid bounds.
 */
private void validate(int row, int col) {
    // Check if the row or column indices are outside the bounds of the grid.
    // The valid range for both 'row' and 'col' is from 0 to gridSize - 1.
    if (row < 0 || row >= gridSize || col < 0 || col >= gridSize) {
        throw new IndexOutOfBoundsException("Index out of bounds.");
    }
}


/**
 * Opens a site at the specified grid coordinates and connects it to its adjacent open neighbors.
 * This method is used to change the state of a site from closed to open. 
 * After opening the site, it will also establish connections to any adjacent open sites
 * to maintain the correct connectivity information in the union-find data structure.
 *
 * @param row The row index of the site to be opened.
 * @param col The column index of the site to be opened.
 */
public void openSite(int row, int col) {
    // Validate the provided indices to ensure they are within the grid bounds.
    validate(row, col);

    // Check if the site at (row, col) is not already open.
    if (!grid[row][col]) {
        // Open the site by setting its value to true.
        grid[row][col] = true;

        // Connect this newly opened site to its adjacent open sites.
        // This is necessary to keep track of full connectivity in the grid.
        // The connections are established in all four directions (up, down, left, right)
        // if those adjacent sites are already open.
        connectToNeighbors(row, col);
    }
}


/**
 * Connects an open site to its adjacent open sites in the grid.
 * This method is used to update the connectivity information in the union-find data structure
 * whenever a new site is opened. It checks the four adjacent sites (up, down, left, right)
 * and if any of them are open, it connects the current site to these open neighbors.
 *
 * @param row The row index of the site to connect.
 * @param col The column index of the site to connect.
 */
private void connectToNeighbors(int row, int col) {
    // Arrays to hold the row and column offsets for adjacent cells.
    // dRow and dCol arrays represent the direction vectors for the adjacent cells.
    int[] dRow = {-1, 1, 0, 0}; // Row offsets for adjacent cells (up, down)
    int[] dCol = {0, 0, -1, 1}; // Column offsets for adjacent cells (left, right)

    // Iterate through each of the four directions (up, down, left, right)
    for (int i = 0; i < dRow.length; i++) {
        // Calculate the coordinates of the adjacent cell in the current direction.
        int newRow = row + dRow[i];
        int newCol = col + dCol[i];

        // Check if the adjacent cell is within the grid bounds and open.
        // If it is open, connect the current site to the adjacent open site.
        // The union operation in the union-find data structure is used to establish this connection.
        // The xyTo1D method converts the 2D coordinates to a 1D index for the union-find structure.
        if (newRow >= 0 && newRow < gridSize && newCol >= 0 && newCol < gridSize && grid[newRow][newCol]) {
            wquFind.union(xyTo1D(row, col), xyTo1D(newRow, newCol));
        }
    }
}


/**
 * Opens sites in the grid randomly based on a given probability.
 * This method iterates through all the sites in the grid and opens each site with a probability 'p'.
 * The decision to open a site is made using a random number generator.
 *
 * @param p The probability with which each site is opened. It's a value between 0 and 1,
 *          where 0 means no site is opened and 1 means all sites are opened.
 */
public void openAllSites(double p) {
    // Iterate through each row of the grid.
    for (int row = 0; row < gridSize; row++) {
        // Iterate through each column of the grid.
        for (int col = 0; col < gridSize; col++) {
            // Generate a random number between 0.0 and 1.0.
            // If the random number is less than the probability 'p', open the site.
            if (StdRandom.uniform() < p) {
                openSite(row, col);
            }
        }
    }
}


/**
 * Checks if the percolation system percolates.
 * A system is said to percolate if there is a path of connected open sites from the top row
 * to the bottom row. This method checks if the virtual top site is connected to the virtual
 * bottom site in the union-find data structure, which would indicate percolation.
 *
 * @return true if the system percolates, false otherwise.
 */
public boolean percolationCheck() {
    // The connected method of the union-find data structure is used to determine if
    // the virtual top site is connected to the virtual bottom site.
    // If they are connected, it means there is a path of open sites from top to bottom,
    // hence the system percolates.
    return wquFind.connected(virtualTop, virtualBottom);
}




/**
 * Displays a grid at a specified position on the canvas, with a white border around it.
 * Each cell of the grid is drawn as a square, with open sites colored differently (e.g., blue)
 * and closed sites in another color (e.g., black). The grid is offset on the canvas based on
 * the provided x and y coordinates.
 *
 * @param xOffset The x-coordinate offset for the grid on the canvas.
 * @param yOffset The y-coordinate offset for the grid on the canvas.
 */
private void displayGrid(int xOffset, int yOffset) {
    int cellSize = 30; // Size of each cell in the grid
    int border = 5; // Width of the white border around the grid

    // Draw a white border around the grid.
    // The filledRectangle method is used to draw the rectangle. The rectangle's center is positioned
    // based on xOffset and yOffset, and its size is calculated to encompass the entire grid plus the border.
    StdDraw.setPenColor(StdDraw.WHITE);
    StdDraw.filledRectangle(xOffset + gridSize * cellSize / 2.0, yOffset + gridSize * cellSize / 2.0, gridSize * cellSize / 2.0 + border, gridSize * cellSize / 2.0 + border);

    // Iterate over each cell of the grid to draw it.
    // Open sites (true in the grid array) are colored blue, and closed sites (false) are colored black.
    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            // Set the pen color based on the site's open or closed status.
            StdDraw.setPenColor(grid[row][col] ? StdDraw.BLUE : StdDraw.BLACK);

            // Draw the cell as a filled square. The position of the square is calculated based on the cell's row and column,
            // adjusted by xOffset and yOffset. The size of each square is half the cellSize, so they fit neatly within the grid.
            StdDraw.filledSquare(xOffset + col * cellSize + cellSize / 2.0, yOffset + (gridSize - row - 1) * cellSize + cellSize / 2.0, cellSize / 2.0);
        }
    }
}


    

    // Main method to run the percolation simulations
    public static void main(String[] args) {
        int numberOfSimulations = 10;
        int n = 10; // Size of each grid
        int cellSize = 30; // Size of each cell in the grid
        int border = 5; // Width of the white border
        int spacing = 10; // Spacing between grids

        // Calculate canvas dimensions for the simulation display
        int gridsPerRow = 5; // Number of grids per row
        int canvasWidth = gridsPerRow * (n * cellSize + 2 * border) + (gridsPerRow - 1) * spacing;
        int numberOfRows = (int) Math.ceil((double) numberOfSimulations / gridsPerRow);
        int canvasHeight = numberOfRows * (n * cellSize + 2 * border) + (numberOfRows - 1) * spacing;

        // Set up the drawing canvas
        StdDraw.setCanvasSize(canvasWidth, canvasHeight);
        StdDraw.setXscale(0, canvasWidth);
        StdDraw.setYscale(0, canvasHeight);
        StdDraw.enableDoubleBuffering();

        // Run simulations and display results
        for (int i = 0; i < numberOfSimulations; i++) {
            double p = StdRandom.uniform(); // Random probability for site opening
            Percolation percolation = new Percolation(n);
            percolation.initializeConnections(); // Establish initial connections
            percolation.openAllSites(p);

            // Calculate the position for each grid
            int xOffset = (i % gridsPerRow) * (n * cellSize + 2 * border + spacing);
            int yOffset = canvasHeight - ((i / gridsPerRow) + 1) * (n * cellSize + 2 * border + spacing);

            // Display each grid and check for percolation
            percolation.displayGrid(xOffset, yOffset);
            boolean doesPercolate = percolation.percolationCheck();
            System.out.println("Simulation " + (i + 1) + ": Probability was " + p + ", the system " + (doesPercolate ? "percolates." : "does not percolate."));
        }

        StdDraw.show(); // Show the final canvas
        


        
    }
}
