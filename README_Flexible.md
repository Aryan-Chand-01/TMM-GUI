# TMM Simulation Tool - Flexible Version

A flexible Transfer Matrix Method (TMM) simulation tool with dynamic layer configuration, dropdown material selection, and automatic material file detection.

## New Features

### Dynamic Layer Configuration
- **Configurable number of layers**: Choose between 1-10 layers using the spinbox
- **Dynamic interface**: Layer input fields are created automatically based on the number selected
- **Flexible structure**: No longer limited to predefined cavity types

### Material Selection
- **Dropdown menus**: Each layer has a dropdown menu to select materials
- **Automatic material detection**: Automatically scans for .txt files and adds them as materials
- **Available materials**:
  - **Metals**: Silver (Ag), Gold (Au)
  - **Dielectrics**: h-BN, PMMA, SiO2, Glass
  - **Active materials**: WS2, MoS2, TDBC
  - **Semiconductors**: Silicon (Si)
  - **Custom materials**: Any .txt file with refractive index data
  - **Others**: Air

### Plot Export Features
- **Copy to clipboard**: Each plot has a "📋 Copy Image" button
- **High-quality export**: Plots are saved as high-DPI PNG images
- **Multiple formats**: Support for PNG, PDF, and SVG export

### User Interface Improvements
- **Scrollable layer panel**: Handle many layers with a scrollable interface
- **Clear labeling**: Each layer is numbered and clearly labeled
- **Real-time updates**: Click "Update Layers" to refresh the interface
- **Material refresh**: Click "🔄 Refresh Materials" to scan for new material files
- **Input validation**: Thickness values are validated as you type

## Adding Custom Materials

You can add new materials by simply placing .txt files in the same directory as the application. The files should have the following format:

```
# Optional comments starting with #
wavelength(µm)  n_real  n_imag
0.400          1.50    0.01
0.450          1.51    0.015
0.500          1.52    0.020
...
```

### Material File Format Requirements:
- **File extension**: Must be .txt
- **Three columns**: wavelength (in µm), real refractive index, imaginary refractive index
- **Numeric data**: All values must be numbers
- **Comments**: Lines starting with # are ignored
- **Filename**: The material name will be the filename without .txt extension

### Example Material Files:
- `SampleMaterial.txt` → Material name: "SampleMaterial"
- `MyCustomGlass.txt` → Material name: "MyCustomGlass"  
- `nSpecialMetal.txt` → Material name: "nSpecialMetal"

After adding new material files:
1. Click the "🔄 Refresh Materials" button
2. The new materials will appear in all dropdown menus
3. A confirmation message will show which materials were found

## How to Use

1. **Set Number of Layers**:
   - Use the spinbox to select the desired number of layers (1-10)
   - Click "Update Layers" to create the input fields

2. **Configure Each Layer**:
   - Select material from the dropdown menu for each layer
   - Enter thickness in nanometers (nm)
   - Layers are numbered from 1 (top) to N (bottom)

3. **Set Substrate**:
   - Choose substrate material from the dropdown

4. **Configure Parameters**:
   - Set wavelength range (min, max, number of points)
   - Set angle range for dispersion calculations
   - Choose polarization (TE, TM, or unpolarized)

5. **Run Calculation**:
   - Click "Calculate" to run the TMM simulation
   - Results will be displayed in plots and analysis panel

## File Structure

- `main_gui.py` - Main GUI with flexible layer system
- `material_properties.py` - Material refractive index calculations
- `tmm_core.py` - TMM calculation engine
- `analysis.py` - Spectral analysis functions
- `plotting.py` - Plotting utilities
- `structure.py` - Structure management (updated for flexible system)
- `run_app.py` - Main application entry point

## Running the Application

```bash
python run_app.py
```

Or use the test script:

```bash
python test_flexible_gui.py
```

## Example Configurations

### Full Microcavity (5 layers):
1. Layer 1: Ag (30 nm) - Bottom mirror
2. Layer 2: hBN (314 nm) - Spacer
3. Layer 3: WS2 (0.7 nm) - Active material  
4. Layer 4: hBN (314 nm) - Spacer
5. Layer 5: Ag (30 nm) - Top mirror

### Simple Structure (3 layers):
1. Layer 1: hBN (100 nm) - Spacer
2. Layer 2: MoS2 (1 nm) - Active material
3. Layer 3: hBN (100 nm) - Spacer

### Single Layer:
1. Layer 1: WS2 (0.7 nm) - Active material only

## Key Improvements

- **Flexibility**: No longer restricted to predefined cavity structures
- **Scalability**: Handle any number of layers from 1 to 10
- **User-friendly**: Intuitive dropdown menus instead of radio buttons
- **Extensible**: Easy to add new materials to the dropdown lists
- **Robust**: Input validation and error handling
- **Scrollable**: Interface scales well with many layers

## Technical Notes

- All thickness values are in nanometers
- Material properties are loaded from the `MaterialProperties` class
- The structure information is displayed in the analysis panel
- Export functionality includes all layer information
- Field calculations account for the flexible layer structure
