"""
Test script for the flexible TMM GUI with dropdown material selection
"""

import tkinter as tk
import sys
import os

# Add the current directory to the path so we can import our modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_gui():
    """Test the flexible GUI"""
    try:
        from main_gui import TMMGUI
        
        print("Starting TMM GUI with flexible layer system...")
        print("Features:")
        print("- Dynamic number of layers (1-10)")
        print("- Dropdown material selection for each layer")
        print("- Flexible thickness input")
        print("- All available materials: Ag, Au, hBN, PMMA, WS2, MoS2, TDBC, SiO2, Silicon, Glass, Air")
        
        root = tk.Tk()
        app = TMMGUI(root)
        
        # Center the window
        root.update_idletasks()
        width = 1400
        height = 900
        x = (root.winfo_screenwidth() // 2) - (width // 2)
        y = (root.winfo_screenheight() // 2) - (height // 2)
        root.geometry(f'{width}x{height}+{x}+{y}')
        
        print("GUI started successfully!")
        print("\nInstructions:")
        print("1. Set the number of layers using the spinbox")
        print("2. Click 'Update Layers' to create layer input fields")
        print("3. Select materials from dropdown menus for each layer")
        print("4. Enter thickness values in nm")
        print("5. Choose substrate material")
        print("6. Set wavelength and angle parameters")
        print("7. Click 'Calculate' to run the simulation")
        
        root.mainloop()
        
    except ImportError as e:
        print(f"Import error: {e}")
        print("Make sure all required modules are available")
    except Exception as e:
        print(f"Error starting GUI: {e}")

if __name__ == "__main__":
    test_gui()
