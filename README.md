# TMM Simulation Tool - Flexible Version

A powerful and user-friendly Transfer Matrix Method (TMM) simulation tool for multilayer thin-film optics, featuring dynamic layer configuration, dropdown material selection, and automatic detection of custom material files.

---

## 🚀 Features

### Dynamic Layer Configuration
- **Configurable number of layers:** Choose between 1–10 layers using a spinbox.
- **Dynamic interface:** Layer input fields are created automatically based on your selection.
- **Flexible structure:** No longer limited to predefined cavity types—model any stack you want!

### Material Selection
- **Dropdown menus:** Each layer has a dropdown menu for material selection.
- **Automatic material detection:** Scans for `.txt` files in the directory and adds them as materials.
- **Available materials:**
  - **Metals:** Silver (Ag), Gold (Au)
  - **Dielectrics:** h-BN, PMMA, SiO₂, Glass
  - **Active materials:** WS₂, MoS₂, TDBC
  - **Semiconductors:** Silicon (Si)
  - **Custom materials:** Any `.txt` file with refractive index data
  - **Others:** Air

### Plot Export and Visualization
- **Copy to clipboard:** Each plot has a 📋 "Copy Image" button.
- **High-quality export:** Save plots as high-DPI PNG, PDF, or SVG.
- **Interactive:** Zoom, pan, and explore results directly in the GUI.

### User Interface Improvements
- **Scrollable layer panel:** Easily handle many layers.
- **Clear labeling:** Each layer is numbered and labeled.
- **Real-time updates:** Click "Update Layers" to refresh the interface.
- **Material refresh:** Click "🔄 Refresh Materials" to scan for new material files.
- **Input validation:** Thickness values are validated as you type.

---

## 🧩 Adding Custom Materials

Add new materials by placing `.txt` files in the application directory. The format should be:

```
# Optional comments starting with #
wavelength(µm)  n_real  n_imag
0.400           1.50    0.01
0.450           1.51    0.015
0.500           1.52    0.020
...
```

**Requirements:**
- File extension: `.txt`
- Three columns: wavelength (µm), real refractive index, imaginary refractive index
- Numeric data only (comments with `#` are ignored)
- Filename (without `.txt`) becomes the material name

**After adding files:**
1. Click "🔄 Refresh Materials"
2. New materials appear in all dropdowns
3. Confirmation message lists detected materials

---

## 🛠 How to Use

1. **Set Number of Layers**
   - Use the spinbox to select 1–10 layers
   - Click "Update Layers" to create input fields

2. **Configure Each Layer**
   - Select material from the dropdown
   - Enter thickness in nanometers (nm)
   - Layers are numbered from 1 (top) to N (bottom)

3. **Set Substrate**
   - Choose substrate material from the dropdown

4. **Configure Parameters**
   - Set wavelength range (min, max, number of points)
   - Set angle range for dispersion calculations (optional)
   - Choose polarization (TE, TM, or unpolarized)

5. **Run Calculation**
   - Click "Calculate" to run the TMM simulation
   - Results are displayed in plots and the analysis panel

---

## 🗂 File Structure

- `main_gui.py` — Main GUI with flexible layer system
- `material_properties.py` — Material refractive index calculations
- `tmm_core.py` — TMM calculation engine
- `analysis.py` — Spectral analysis functions
- `plotting.py` — Plotting utilities
- `structure.py` — Structure management (flexible system)
- `run_app.py` — Main application entry point

---

## ▶️ Running the Application

```bash
python run_app.py
```

Or use the test script:

```bash
python test_flexible_gui.py
```

---


## 🌟 Key Improvements

- **Flexibility:** No restriction to predefined cavity structures
- **Scalability:** Model any stack from 1 to 10 layers
- **User-friendly:** Intuitive dropdowns and real-time feedback
- **Extensible:** Add new materials with a simple file drop
- **Robust:** Input validation and error handling
- **Scrollable:** Interface scales for many layers

---

## ⚙️ Technical Notes

- All thickness values are in nanometers (nm)
- Material properties are loaded from the `MaterialProperties` class
- Structure information is displayed in the analysis panel
- Export functionality includes all layer information
- Field calculations account for the flexible layer structure

---

## 📚 License and Citation

This tool is provided for academic and research use. If you use it in your work, please cite appropriately.

---

## 👨‍💻 Author

Developed by Aryan Chand, Indian Institute of Space Science and Technology, Trivandrum.

---

## 💡 Contributing


