# Map Projections with Gnuplot

## Environment
- Gnuplot 5.2 or newer (arrays are used)

## Setup
1. Clone the repository and `cd` into the repository root (the same directory as `world.plt`).
2. Download `world_110m.txt` from the gnuplotting site and place it next to `world.plt` (e.g. using `curl`).
   ```sh
   curl -O https://example.com/path/to/world_110m.txt
   ```
3. Start Gnuplot and call the script:
   - Interactive:
     ```gnuplot
     call 'world.plt' makegraticule
     ```
   - Or run in batch mode:
     ```sh
     gnuplot -c world.plt makegraticule
     ```

## Usage
1. In Gnuplot, load the script to display the world map:
   ```gnuplot
   load 'world.plt'
   ```
   The map will be drawn.
2. To see the commands that generated the map, run:
   ```gnuplot
   show plot
   ```
3. The available projection functions are defined inside `world.plt`; inspect the file for the list.

## Distortion evaluation
1. Load `world.plt` in Gnuplot as above.
2. Define a projection function with a single argument (fix any parameters).
3. Evaluate distortion:
   ```gnuplot
   call 'world.plt' defstrainfunc <projection_function>
   call 'world.plt' evaluate <projection_function> [<metric_function>]
   ```
4. Note that `makegraticule` only creates a mesh covering the northern‑eastern quarter of the Earth. For asymmetric projections you may need to rotate the projection to evaluate the southern hemisphere and average the results.
5. To test convergence with respect to grid resolution:
   ```gnuplot
   call 'world.plt' defstrainfunc <projection_function>
   call 'world.plt' gridsizeconvergence <metric_function>
   ```
   This will output results for grid spacings from 2° down to 0.2°. Some projections may take several minutes.

## Platform notes
- The instructions work on both Linux and Windows. On Windows ensure that the Gnuplot executable is in your `PATH` or use the full path to `gnuplot.exe`.
