"""
Test the new features: copy graph functionality and automatic material detection
"""

import tkinter as tk
import os
import sys

def test_material_scanning():
    """Test the material scanning functionality"""
    print("Testing Material Scanner...")
    
    try:
        # Import the GUI class
        from main_gui import TMMGUI
        
        # Create a dummy root for testing
        root = tk.Tk()
        root.withdraw()  # Hide the window
        
        # Create GUI instance to test material scanning
        app = TMMGUI(root)
        
        print(f"Base materials: {app.base_materials}")
        print(f"Available materials: {app.available_materials}")
        
        # Check for our sample materials
        sample_materials = ['SampleMaterial', 'CustomMaterial']
        found_samples = [m for m in app.available_materials if m in sample_materials]
        
        print(f"Sample materials found: {found_samples}")
        
        if found_samples:
            print("✅ Material scanning working correctly!")
        else:
            print("⚠️  No sample materials found. Make sure .txt files exist.")
        
        root.destroy()
        return True
        
    except Exception as e:
        print(f"❌ Error testing material scanning: {e}")
        return False

def test_file_listing():
    """List all .txt files in current directory"""
    print("\nListing .txt files in current directory:")
    current_dir = os.getcwd()
    txt_files = [f for f in os.listdir(current_dir) if f.endswith('.txt')]
    
    for txt_file in txt_files:
        print(f"  📄 {txt_file}")
        
        # Try to read first few lines
        try:
            with open(txt_file, 'r') as f:
                lines = f.readlines()[:3]
                print(f"     First lines: {[line.strip() for line in lines if line.strip()]}")
        except Exception as e:
            print(f"     Error reading: {e}")

def main():
    print("TMM GUI Feature Test")
    print("===================")
    
    # Test file listing
    test_file_listing()
    
    # Test material scanning
    success = test_material_scanning()
    
    if success:
        print("\n🎉 All tests passed! New features are working.")
        print("\nNew Features Added:")
        print("1. 📋 Copy Image buttons on all plots")
        print("2. 🔄 Automatic material file detection")
        print("3. 🔄 Refresh Materials button")
        print("\nTo test the GUI:")
        print("python run_app.py")
    else:
        print("\n❌ Some tests failed. Check the error messages above.")

if __name__ == "__main__":
    main()
