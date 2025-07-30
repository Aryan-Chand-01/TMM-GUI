"""
TMM Simulation Tool - Entry Point (Flexible Version)
Run this file to start the application with flexible layer configuration
"""

import tkinter as tk
from main_gui import TMMGUI

def main():
    """Main entry point for the TMM simulation tool"""
    print("TMM Simulation Tool - Flexible Version with New Features")
    print("=" * 60)
    print("Features:")
    print("- Dynamic layer configuration (1-10 layers)")
    print("- Dropdown material selection")
    print("- Flexible thickness input")
    print("- Built-in materials: Ag, Au, hBN, PMMA, WS2, MoS2, TDBC, SiO2, Silicon, Glass, Air")
    print("- 📋 Copy Image buttons on all plots")
    print("- 🔄 Automatic detection of .txt material files")
    print("- 🔄 Refresh Materials button")
    print("\nNew .txt files in the directory will automatically be detected as materials!")
    print("File format: wavelength(µm) n_real n_imag")
    print("\nStarting GUI...")
    
    root = tk.Tk()
    app = TMMGUI(root)
    
    # Print detected materials
    print(f"\nDetected {len(app.available_materials)} materials:")
    for i, material in enumerate(app.available_materials, 1):
        print(f"  {i:2d}. {material}")
    
    # Center the window on screen
    root.update_idletasks()
    width = root.winfo_width()
    height = root.winfo_height()
    x = (root.winfo_screenwidth() // 2) - (width // 2)
    y = (root.winfo_screenheight() // 2) - (height // 2)
    root.geometry(f'{width}x{height}+{x}+{y}')
    
    print("\n🎉 GUI started! Look for:")
    print("- 📋 Copy Image buttons next to plot toolbars")
    print("- 🔄 Refresh Materials button in Structure tab")
    
    root.mainloop()

if __name__ == "__main__":
    main()