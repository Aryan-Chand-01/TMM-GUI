import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import threading
import os
import glob
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import io

try:
    from PIL import Image
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False

# Import our modules
from tmm_core import TMMCore
from analysis import TMMAnalysis
from plotting import TMMPlotter
from structure import StructureManager

class TMMGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("TMM Simulation Tool - Flexible Version")
        self.root.geometry("1400x900")
        
        # Initialize components
        self.tmm = TMMCore()
        self.analyzer = TMMAnalysis()
        self.plotter = TMMPlotter()
        self.structure_manager = StructureManager()
        self.current_results = None
        
        # Available materials (base materials + scanned files)
        self.base_materials = [
            "Air", "Ag", "Au", "hBN", "PMMA", "WS2", "MoS2", "TDBC", 
            "SiO2", "Silicon", "Glass"
        ]
        self.available_materials = self.scan_for_material_files()
        
        # Layer configuration
        self.num_layers = tk.IntVar(value=3)
        self.layer_materials = []
        self.layer_thicknesses = []
        self.layer_frame = None
        
        # Create main frames
        self.create_widgets()
        self.setup_default_values()
    
    def scan_for_material_files(self):
        """Scan current directory for .txt files and add them as materials"""
        materials = self.base_materials.copy()
        
        try:
            # Get current directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            
            # Look for .txt files (excluding specific system files)
            exclude_files = ['README.txt', 'LICENSE.txt', 'requirements.txt', 'config.txt']
            
            for txt_file in glob.glob(os.path.join(current_dir, "*.txt")):
                filename = os.path.basename(txt_file)
                
                # Skip excluded files
                if filename.lower() in [f.lower() for f in exclude_files]:
                    continue
                
                # Extract material name (remove .txt extension)
                material_name = os.path.splitext(filename)[0]
                
                # Check if file has the expected format (at least 3 columns of numbers)
                try:
                    with open(txt_file, 'r') as f:
                        lines = f.readlines()
                        # Check first few lines for numeric data
                        valid_file = False
                        for line in lines[:5]:  # Check first 5 lines
                            parts = line.strip().split()
                            if len(parts) >= 3:
                                try:
                                    float(parts[0])  # wavelength
                                    float(parts[1])  # n_real
                                    float(parts[2])  # n_imag
                                    valid_file = True
                                    break
                                except ValueError:
                                    continue
                        
                        if valid_file and material_name not in self.base_materials:
                            materials.append(material_name)
                            print(f"Found material file: {filename} -> Added '{material_name}' to materials list")
                
                except Exception as e:
                    print(f"Error reading {filename}: {e}")
            
        except Exception as e:
            print(f"Error scanning for material files: {e}")
        
        return materials
    
    def copy_plot_to_clipboard(self, figure):
        """Copy plot to clipboard as image"""
        try:
            if PIL_AVAILABLE:
                # Save figure to BytesIO buffer
                buf = io.BytesIO()
                figure.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                buf.seek(0)
                
                # Convert to PIL Image
                img = Image.open(buf)
                
                # Copy to clipboard
                self.root.clipboard_clear()
                
                # For Windows, we need to save to a temporary file and copy
                import tempfile
                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as temp_file:
                    img.save(temp_file.name)
                    temp_path = temp_file.name
                
                # Copy file path to clipboard (user can paste in any application)
                self.root.clipboard_append(temp_path)
                messagebox.showinfo("Copy Complete", f"Plot image saved to temporary file.\nPath copied to clipboard: {temp_path}")
                
            else:
                # Fallback: save to file
                filename = filedialog.asksaveasfilename(
                    defaultextension=".png",
                    filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("SVG files", "*.svg")]
                )
                if filename:
                    figure.savefig(filename, dpi=300, bbox_inches='tight')
                    messagebox.showinfo("Save Complete", f"Plot saved to {filename}")
                    
        except Exception as e:
            messagebox.showerror("Copy Error", f"Error copying plot: {str(e)}")
        
    def create_widgets(self):
        # Main container
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Left panel for controls
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        
        # Right panel for plots
        plot_frame = ttk.Frame(main_frame)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        self.create_control_panel(control_frame)
        self.create_plot_panel(plot_frame)
    
    def create_control_panel(self, parent):
        # Control panel with notebook for organization
        notebook = ttk.Notebook(parent)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Structure tab
        structure_frame = ttk.Frame(notebook)
        notebook.add(structure_frame, text="Structure")
        
        # Parameters tab
        params_frame = ttk.Frame(notebook)
        notebook.add(params_frame, text="Parameters")
        
        # Analysis tab
        analysis_frame = ttk.Frame(notebook)
        notebook.add(analysis_frame, text="Analysis")
        
        self.create_structure_controls(structure_frame)
        self.create_parameter_controls(params_frame)
        self.create_analysis_controls(analysis_frame)
    
    def create_structure_controls(self, parent):
        # Number of layers selection
        layers_config_frame = ttk.LabelFrame(parent, text="Layer Configuration")
        layers_config_frame.pack(fill=tk.X, padx=5, pady=5)
        
        config_inner = ttk.Frame(layers_config_frame)
        config_inner.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(config_inner, text="Number of Layers:").grid(row=0, column=0, sticky=tk.W, pady=2)
        num_layers_spinbox = ttk.Spinbox(config_inner, from_=1, to=10, width=10, 
                                        textvariable=self.num_layers,
                                        command=self.update_layer_inputs)
        num_layers_spinbox.grid(row=0, column=1, padx=5, pady=2)
        
        # Button to apply changes
        ttk.Button(config_inner, text="Update Layers", 
                  command=self.update_layer_inputs).grid(row=0, column=2, padx=10, pady=2)
        
        # Refresh materials button
        ttk.Button(config_inner, text="🔄 Refresh Materials", 
                  command=self.refresh_materials).grid(row=0, column=3, padx=5, pady=2)
        
        # Substrate selection
        substrate_frame = ttk.LabelFrame(parent, text="Substrate")
        substrate_frame.pack(fill=tk.X, padx=5, pady=5)
        
        substrate_inner = ttk.Frame(substrate_frame)
        substrate_inner.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(substrate_inner, text="Material:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.substrate = tk.StringVar(value="Air")
        self.substrate_combo = ttk.Combobox(substrate_inner, textvariable=self.substrate, 
                                           values=self.available_materials, 
                                           state="readonly", width=15)
        self.substrate_combo.grid(row=0, column=1, padx=5, pady=2)
        
        # Container for dynamic layer inputs
        self.layer_container = ttk.Frame(parent)
        self.layer_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Initialize with default layers
        self.update_layer_inputs()
    
    def refresh_materials(self):
        """Refresh the list of available materials by scanning for new txt files"""
        old_materials = self.available_materials.copy()
        self.available_materials = self.scan_for_material_files()
        
        # Update substrate combobox
        self.substrate_combo['values'] = self.available_materials
        
        # Update layer material comboboxes if they exist
        if hasattr(self, 'layer_material_combos'):
            for combo in self.layer_material_combos:
                combo['values'] = self.available_materials
        
        # Show message about new materials found
        new_materials = [m for m in self.available_materials if m not in old_materials]
        if new_materials:
            messagebox.showinfo("Materials Refreshed", 
                              f"Found {len(new_materials)} new material(s):\n" + 
                              "\n".join(new_materials))
        else:
            messagebox.showinfo("Materials Refreshed", "No new materials found.")
    
    def update_layer_inputs(self):
        """Update the layer input fields based on number of layers"""
        # Clear existing layer frame
        if self.layer_frame:
            self.layer_frame.destroy()
        
        # Clear existing variables
        self.layer_materials.clear()
        self.layer_thicknesses.clear()
        
        # Initialize list to store material comboboxes for refresh functionality
        self.layer_material_combos = []
        
        # Create new frame for layers
        self.layer_frame = ttk.LabelFrame(self.layer_container, text="Layer Properties")
        self.layer_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create scrollable frame for many layers
        canvas = tk.Canvas(self.layer_frame, height=300)
        scrollbar = ttk.Scrollbar(self.layer_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Header
        header_frame = ttk.Frame(scrollable_frame)
        header_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(header_frame, text="Layer #", width=8, font=("Arial", 9, "bold")).grid(row=0, column=0, padx=2)
        ttk.Label(header_frame, text="Material", width=12, font=("Arial", 9, "bold")).grid(row=0, column=1, padx=2)
        ttk.Label(header_frame, text="Thickness (nm)", width=15, font=("Arial", 9, "bold")).grid(row=0, column=2, padx=2)
        
        # Create input rows for each layer
        num_layers = self.num_layers.get()
        default_materials = ["Ag", "hBN", "WS2", "hBN", "Ag"]
        default_thicknesses = [30.0, 314.0, 0.7, 314.0, 30.0]
        
        for i in range(num_layers):
            layer_frame = ttk.Frame(scrollable_frame)
            layer_frame.pack(fill=tk.X, padx=5, pady=2)
            
            # Layer number
            ttk.Label(layer_frame, text=f"{i+1}", width=8).grid(row=0, column=0, padx=2)
            
            # Material dropdown
            material_var = tk.StringVar()
            if i < len(default_materials):
                material_var.set(default_materials[i])
            else:
                material_var.set("Air")
            
            material_combo = ttk.Combobox(layer_frame, textvariable=material_var,
                                        values=self.available_materials,
                                        state="readonly", width=12)
            material_combo.grid(row=0, column=1, padx=2)
            self.layer_materials.append(material_var)
            self.layer_material_combos.append(material_combo)  # Store reference for refresh
            
            # Thickness entry
            thickness_var = tk.DoubleVar()
            if i < len(default_thicknesses):
                thickness_var.set(default_thicknesses[i])
            else:
                thickness_var.set(10.0)
            
            thickness_entry = ttk.Entry(layer_frame, textvariable=thickness_var, width=15)
            thickness_entry.grid(row=0, column=2, padx=2)
            self.layer_thicknesses.append(thickness_var)
            
            # Add validation
            def validate_number(value):
                try:
                    float(value)
                    return True
                except ValueError:
                    return False
            
            vcmd = (self.root.register(validate_number), '%P')
            thickness_entry.config(validate='key', validatecommand=vcmd)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Bind mousewheel to canvas
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        

    def create_parameter_controls(self, parent):
        # Wavelength range
        wl_frame = ttk.LabelFrame(parent, text="Wavelength Range")
        wl_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(wl_frame, text="Min (nm):").grid(row=0, column=0)
        self.wl_min = tk.DoubleVar(value=476)
        ttk.Entry(wl_frame, textvariable=self.wl_min, width=10).grid(row=0, column=1)
        
        ttk.Label(wl_frame, text="Max (nm):").grid(row=1, column=0)
        self.wl_max = tk.DoubleVar(value=714)
        ttk.Entry(wl_frame, textvariable=self.wl_max, width=10).grid(row=1, column=1)
        
        ttk.Label(wl_frame, text="Points:").grid(row=2, column=0)
        self.wl_points = tk.IntVar(value=200)
        ttk.Entry(wl_frame, textvariable=self.wl_points, width=10).grid(row=2, column=1)
        
        # Angle range
        angle_frame = ttk.LabelFrame(parent, text="Angle Range")
        angle_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(angle_frame, text="Min (deg):").grid(row=0, column=0)
        self.angle_min = tk.DoubleVar(value=0)
        ttk.Entry(angle_frame, textvariable=self.angle_min, width=10).grid(row=0, column=1)
        
        ttk.Label(angle_frame, text="Max (deg):").grid(row=1, column=0)
        self.angle_max = tk.DoubleVar(value=85)
        ttk.Entry(angle_frame, textvariable=self.angle_max, width=10).grid(row=1, column=1)
        
        ttk.Label(angle_frame, text="Step:").grid(row=2, column=0)
        self.angle_step = tk.DoubleVar(value=2)
        ttk.Entry(angle_frame, textvariable=self.angle_step, width=10).grid(row=2, column=1)
        
        # Polarization
        pol_frame = ttk.LabelFrame(parent, text="Polarization")
        pol_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.polarization = tk.StringVar(value="")
        ttk.Radiobutton(pol_frame, text="Unpolarized", variable=self.polarization, 
                       value="").pack(anchor=tk.W)
        ttk.Radiobutton(pol_frame, text="TE", variable=self.polarization, 
                       value="TE").pack(anchor=tk.W)
        ttk.Radiobutton(pol_frame, text="TM", variable=self.polarization, 
                       value="TM").pack(anchor=tk.W)
        
        # Options
        options_frame = ttk.LabelFrame(parent, text="Options")
        options_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.calc_field = tk.BooleanVar(value=True)
        ttk.Checkbutton(options_frame, text="Calculate Electric Field", 
                       variable=self.calc_field).pack(anchor=tk.W)
        
        self.calc_dispersion = tk.BooleanVar(value=False)
        ttk.Checkbutton(options_frame, text="Calculate Dispersion", 
                       variable=self.calc_dispersion).pack(anchor=tk.W)
        
        # Calculate button
        calculate_btn = ttk.Button(parent, text="Calculate", command=self.run_calculation)
        calculate_btn.pack(pady=10)
    
    def create_analysis_controls(self, parent):
        # Results display
        results_frame = ttk.LabelFrame(parent, text="Results")
        results_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.results_text = tk.Text(results_frame, height=15, width=40)
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scrollbar.set)
        
        self.results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Export button
        ttk.Button(parent, text="Export Data", command=self.export_data).pack(pady=5)
    
    def create_plot_panel(self, parent):
        # Create notebook for different plot types
        self.plot_notebook = ttk.Notebook(parent)
        self.plot_notebook.pack(fill=tk.BOTH, expand=True)
        
        # Transmission/Reflection plot
        self.tr_frame = ttk.Frame(self.plot_notebook)
        self.plot_notebook.add(self.tr_frame, text="T/R/A")
        
        # Field plot
        self.field_frame = ttk.Frame(self.plot_notebook)
        self.plot_notebook.add(self.field_frame, text="Field")
        
        # Combined dispersion plot (renamed from "All Dispersions")
        self.disp_multi_frame = ttk.Frame(self.plot_notebook)
        self.plot_notebook.add(self.disp_multi_frame, text="Dispersion Plot")
        
        # Create matplotlib figures
        self.create_plot_figures()
    
    def create_plot_figures(self):
        # T/R/A plot
        self.tr_fig = Figure(figsize=(8, 6), dpi=100)
        self.tr_ax = self.tr_fig.add_subplot(111)
        self.tr_canvas = FigureCanvasTkAgg(self.tr_fig, self.tr_frame)
        self.tr_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar and copy button for T/R/A plot
        toolbar_frame1 = tk.Frame(self.tr_frame)
        toolbar_frame1.pack(side=tk.BOTTOM, fill=tk.X)
        self.tr_toolbar = NavigationToolbar2Tk(self.tr_canvas, toolbar_frame1)
        self.tr_toolbar.update()
        
        # Add checkboxes for T/R/A selection
        self.show_trans_var = tk.BooleanVar(value=True)
        self.show_refl_var = tk.BooleanVar(value=True)
        self.show_abs_var = tk.BooleanVar(value=True)
        cb_frame = tk.Frame(toolbar_frame1)
        cb_frame.pack(side=tk.LEFT, padx=5)
        ttk.Checkbutton(cb_frame, text="Transmission", variable=self.show_trans_var, command=self.update_tra_plot).pack(side=tk.LEFT)
        ttk.Checkbutton(cb_frame, text="Reflection", variable=self.show_refl_var, command=self.update_tra_plot).pack(side=tk.LEFT)
        ttk.Checkbutton(cb_frame, text="Absorption", variable=self.show_abs_var, command=self.update_tra_plot).pack(side=tk.LEFT)

        # Add copy button for T/R/A plot
        copy_btn1 = ttk.Button(toolbar_frame1, text="📋 Copy Image", 
                              command=lambda: self.copy_plot_to_clipboard(self.tr_fig))
        copy_btn1.pack(side=tk.RIGHT, padx=5)
        
        # Field plot
        self.field_fig = Figure(figsize=(8, 6), dpi=100)
        self.field_ax = self.field_fig.add_subplot(111)
        self.field_canvas = FigureCanvasTkAgg(self.field_fig, self.field_frame)
        self.field_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar and copy button for field plot
        toolbar_frame2 = tk.Frame(self.field_frame)
        toolbar_frame2.pack(side=tk.BOTTOM, fill=tk.X)
        self.field_toolbar = NavigationToolbar2Tk(self.field_canvas, toolbar_frame2)
        self.field_toolbar.update()
        
        # Add copy button for field plot
        copy_btn2 = ttk.Button(toolbar_frame2, text="📋 Copy Image", 
                              command=lambda: self.copy_plot_to_clipboard(self.field_fig))
        copy_btn2.pack(side=tk.RIGHT, padx=5)
        
        # Combined dispersion plot
        self.disp_multi_fig = Figure(figsize=(12, 10), dpi=100)
        self.disp_multi_canvas = FigureCanvasTkAgg(self.disp_multi_fig, self.disp_multi_frame)
        self.disp_multi_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add toolbar and copy button for dispersion plot
        toolbar_frame3 = tk.Frame(self.disp_multi_frame)
        toolbar_frame3.pack(side=tk.BOTTOM, fill=tk.X)
        self.disp_multi_toolbar = NavigationToolbar2Tk(self.disp_multi_canvas, toolbar_frame3)
        self.disp_multi_toolbar.update()
        
        # Add copy button for dispersion plot
        copy_btn3 = ttk.Button(toolbar_frame3, text="📋 Copy Image", 
                              command=lambda: self.copy_plot_to_clipboard(self.disp_multi_fig))
        copy_btn3.pack(side=tk.RIGHT, padx=5)
    
    def setup_default_values(self):
        pass
    
    def get_structure(self):
        """Get current structure configuration from dynamic layers"""
        try:
            materials = []
            thicknesses = []
            
            # Get data from dynamic layer inputs
            for i in range(len(self.layer_materials)):
                material = self.layer_materials[i].get()
                thickness = self.layer_thicknesses[i].get()
                
                materials.append(material)
                thicknesses.append(thickness)
            
            return thicknesses, materials
        
        except Exception as e:
            print(f"Error getting structure: {e}")
            # Return default structure if error
            return [30.0, 314.0, 0.7], ["Ag", "hBN", "WS2"]

    def clear_all_plots(self):
        """Clear all plots before new calculation"""
        try:
            # Clear all stored colorbars in plotter first
            self.plotter.clear_all_colorbars()
            
            # Clear T/R/A plot
            self.tr_ax.clear()
            self.tr_ax.text(0.5, 0.5, 'Calculating...', 
                           transform=self.tr_ax.transAxes, ha='center', va='center', fontsize=12)
            self.tr_canvas.draw()
            
            # Clear field plot completely (including colorbars)
            self.field_fig.clear()
            self.field_ax = self.field_fig.add_subplot(111)
            self.field_ax.text(0.5, 0.5, 'Calculating...', 
                              transform=self.field_ax.transAxes, ha='center', va='center', fontsize=12)
            self.field_canvas.draw()
            
            # Clear dispersion plot completely (including colorbars)
            self.disp_multi_fig.clear()
            ax = self.disp_multi_fig.add_subplot(111)
            ax.text(0.5, 0.5, 'Calculating...', 
                   transform=ax.transAxes, ha='center', va='center', fontsize=12)
            self.disp_multi_canvas.draw()
            
        except Exception as e:
            print(f"Error clearing plots: {e}")

    
    def run_calculation(self):
        """Run TMM calculation in separate thread"""
        # Clear all plots immediately when Calculate is clicked
        self.clear_all_plots()
        
        # Clear results text
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(1.0, "Calculating...\nPlease wait...")
        
        def calculate():
            try:
                # Get parameters
                wavelengths = np.linspace(self.wl_min.get(), self.wl_max.get(), self.wl_points.get())
                
                if self.calc_dispersion.get():
                    angles = np.arange(self.angle_min.get(), 
                                     self.angle_max.get() + self.angle_step.get(), 
                                     self.angle_step.get())
                else:
                    angles = np.array([0])  # Normal incidence only
                
                thicknesses, materials = self.get_structure()
                substrate = self.substrate.get()
                polarization = self.polarization.get()
                calc_field = self.calc_field.get()
                
                # Run calculation
                self.current_results = self.tmm.tmm_calculation(
                    wavelengths, thicknesses, materials, substrate, 
                    polarization, angles, calc_field
                )
                
                # Store additional info
                self.current_results['wavelengths'] = wavelengths
                self.current_results['angles'] = angles
                self.current_results['thicknesses'] = thicknesses
                self.current_results['materials'] = materials
                
                # Update plots on main thread
                self.root.after(0, self.update_plots)
                
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Calculation Error", str(e)))
                # Clear the "Calculating..." message on error
                self.root.after(0, self.clear_calculating_messages)
        
        # Run in thread to keep GUI responsive
        thread = threading.Thread(target=calculate)
        thread.daemon = True
        thread.start()
    
    def clear_calculating_messages(self):
        """Clear 'Calculating...' messages from plots"""
        try:
            # Clear all axes and show empty plots
            self.tr_ax.clear()
            self.tr_ax.set_title('Transmission, Reflection, and Absorption')
            self.tr_ax.set_xlabel('Wavelength (nm)')
            self.tr_ax.set_ylabel('T, R, A')
            self.tr_canvas.draw()
            
            self.field_ax.clear()
            self.field_ax.set_title('Electric Field Intensity')
            self.field_ax.set_xlabel('Wavelength (nm)')
            self.field_ax.set_ylabel('Position (nm)')
            self.field_canvas.draw()
            
            # Clear dispersion plot
            self.disp_multi_fig.clear()
            self.disp_multi_canvas.draw()
            
        except Exception as e:
            print(f"Error clearing calculating messages: {e}")
    
    def update_plots(self):
        """Update all plots with current results"""
        if self.current_results is None:
            return
        
        # Update T/R/A plot
        self.update_tra_plot()
        
        # Update field plot if available
        if 'field' in self.current_results:
            self.update_field_plot()
        
        # Update dispersion plot if calculated
        if len(self.current_results['angles']) > 1:
            self.update_dispersion_plot()
        
        # Update results text
        self.update_results_text()
    
    def update_tra_plot(self):
        """Update transmission/reflection/absorption plot"""
        wavelengths = self.current_results['wavelengths']
        transmission = self.current_results['transmission'][:, 0]
        reflection = self.current_results['reflection'][:, 0]
        absorption = self.current_results['absorption'][:, 0]

        show_trans = self.show_trans_var.get() if hasattr(self, 'show_trans_var') else True
        show_refl = self.show_refl_var.get() if hasattr(self, 'show_refl_var') else True
        show_abs = self.show_abs_var.get() if hasattr(self, 'show_abs_var') else True

        self.plotter.plot_tra_spectrum(
            self.tr_fig, self.tr_ax, wavelengths,
            transmission, reflection, absorption,
            show_trans, show_refl, show_abs
        )
        self.tr_canvas.draw()
    
    def update_field_plot(self):
        """Update electric field plot"""
        if 'field' not in self.current_results:
            return
        
        wavelengths = self.current_results['wavelengths']
        z = self.current_results['z'] * 1e9  # Convert to nm
        field = self.current_results['field']
        thicknesses = self.current_results['thicknesses']
        
        self.plotter.plot_field_intensity(
            self.field_fig, self.field_ax, wavelengths, 
            z, field, thicknesses
        )
        self.field_canvas.draw()
    
    def update_dispersion_plot(self):
        """Update the combined dispersion plot"""
        if len(self.current_results['angles']) <= 1:
            # Clear dispersion plot if no angle sweep
            self.disp_multi_fig.clear()
            ax = self.disp_multi_fig.add_subplot(111)
            ax.text(0.5, 0.5, 'Enable "Calculate Dispersion"\nto see angle-resolved plots', 
                   transform=ax.transAxes, ha='center', va='center', fontsize=12)
            self.disp_multi_canvas.draw()
            return
        
        try:
            wavelengths = self.current_results['wavelengths']
            angles = self.current_results['angles']
            transmission = self.current_results['transmission']
            reflection = self.current_results['reflection']
            absorption = self.current_results['absorption']
            
            # Update combined dispersion plot
            self.plotter.plot_multiple_dispersion(
                self.disp_multi_fig, wavelengths, angles, 
                transmission, reflection, absorption
            )
            self.disp_multi_canvas.draw()
            
        except Exception as e:
            print(f"Error updating dispersion plot: {e}")
            # Show error in plot
            self.disp_multi_fig.clear()
            ax = self.disp_multi_fig.add_subplot(111)
            ax.text(0.5, 0.5, f'Error plotting dispersion:\n{str(e)}', 
                   transform=ax.transAxes, ha='center', va='center')
            self.disp_multi_canvas.draw()
    
    def update_results_text(self):
        """Update results text display"""
        self.results_text.delete(1.0, tk.END)
        
        if self.current_results is None:
            self.results_text.insert(1.0, "No calculation results available.\nRun calculation to see analysis.")
            return
        
        try:
            # Get structure info using new flexible system
            materials = self.current_results['materials']
            thicknesses = self.current_results['thicknesses']
            
            text = self.structure_manager.get_flexible_structure_info(materials, thicknesses)
            text += f"Substrate: {self.substrate.get()}\n"
            text += f"Polarization: {self.polarization.get() if self.polarization.get() else 'Unpolarized'}\n\n"
            
            # Analysis
            wavelengths = self.current_results['wavelengths']
            transmission = self.current_results['transmission']
            reflection = self.current_results['reflection']
            absorption = self.current_results['absorption']
            
            # Use first angle data for analysis
            trans_1d = transmission[:, 0] if transmission.ndim > 1 else transmission
            refl_1d = reflection[:, 0] if reflection.ndim > 1 else reflection
            abs_1d = absorption[:, 0] if absorption.ndim > 1 else absorption
            
            analysis = self.analyzer.analyze_spectrum(wavelengths, trans_1d, refl_1d, abs_1d)
            
            text += "Analysis Results:\n"
            text += f"Max Transmission: {analysis['max_transmission']:.4f}\n"
            text += f"Min Reflection: {analysis['min_reflection']:.4f}\n"
            text += f"Max Absorption: {analysis['max_absorption']:.4f}\n\n"
            
            # Wavelength at max transmission
            max_trans_idx = np.argmax(trans_1d)
            text += f"Max Transmission at: {wavelengths[max_trans_idx]:.1f} nm\n"
            
            # Wavelength at min reflection
            min_refl_idx = np.argmin(refl_1d)
            text += f"Min Reflection at: {wavelengths[min_refl_idx]:.1f} nm\n\n"
            
            # Rabi splitting
            if analysis['rabi_splitting']['found']:
                rs = analysis['rabi_splitting']
                text += "Rabi Splitting Found:\n"
                text += f"Peak 1: {rs['peak1_nm']:.1f} nm ({rs['peak1_meV']:.1f} meV)\n"
                text += f"Peak 2: {rs['peak2_nm']:.1f} nm ({rs['peak2_meV']:.1f} meV)\n"
                text += f"Splitting: {rs['splitting_nm']:.1f} nm ({rs['splitting_meV']:.1f} meV)\n\n"
            else:
                text += "No clear Rabi splitting detected\n\n"
            
            # FWHM
            if analysis['fwhm']['found']:
                fwhm = analysis['fwhm']
                text += f"FWHM: {fwhm['fwhm_nm']:.1f} nm ({fwhm['fwhm_meV']:.1f} meV)\n"
                text += f"Q-factor: {fwhm['q_factor']:.1f}\n"
                text += f"Center: {fwhm['center_nm']:.1f} nm\n\n"
            else:
                text += "FWHM could not be determined\n\n"
            
            # Peak/dip information
            peaks_dips = analysis['peaks_dips']
            if len(peaks_dips['peaks']) > 0:
                text += f"Transmission peaks found at: {[f'{p:.1f}' for p in peaks_dips['peaks']]} nm\n"
            if len(peaks_dips['dips']) > 0:
                text += f"Reflection dips found at: {[f'{d:.1f}' for d in peaks_dips['dips']]} nm\n"
            
            # Dispersion info
            if len(self.current_results['angles']) > 1:
                text += f"\nDispersion calculated over {len(self.current_results['angles'])} angles\n"
                text += f"Angle range: {self.current_results['angles'][0]:.1f}° to {self.current_results['angles'][-1]:.1f}°\n"
            
            self.results_text.insert(1.0, text)
            
        except Exception as e:
            error_text = f"Error updating analysis:\n{str(e)}\n\nRaw data available:\n"
            if self.current_results:
                for key in self.current_results.keys():
                    if hasattr(self.current_results[key], 'shape'):
                        error_text += f"{key}: {self.current_results[key].shape}\n"
                    else:
                        error_text += f"{key}: {type(self.current_results[key])}\n"
            self.results_text.insert(1.0, error_text)
    
    def export_data(self):
        """Export calculation results to CSV"""
        if self.current_results is None:
            messagebox.showwarning("No Data", "No calculation results to export")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                import csv
                
                # Create data for CSV export
                wavelengths = self.current_results['wavelengths']
                transmission = self.current_results['transmission'][:, 0]
                reflection = self.current_results['reflection'][:, 0]
                absorption = self.current_results['absorption'][:, 0]
                
                # Write CSV file
                with open(filename, 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    
                    # Prepare header and data
                    header = ['Wavelength_nm', 'Transmission', 'Reflection', 'Absorption']
                    
                    # If dispersion data exists, add angle columns
                    if len(self.current_results['angles']) > 1:
                        angles = self.current_results['angles']
                        for angle in angles:
                            header.extend([f'Transmission_{angle:.1f}deg', f'Reflection_{angle:.1f}deg', f'Absorption_{angle:.1f}deg'])
                    
                    writer.writerow(header)
                    
                    # Write data rows
                    for i in range(len(wavelengths)):
                        row = [wavelengths[i], transmission[i], reflection[i], absorption[i]]
                        
                        # Add angle-resolved data if available
                        if len(self.current_results['angles']) > 1:
                            for j, angle in enumerate(self.current_results['angles']):
                                row.extend([
                                    self.current_results['transmission'][i, j],
                                    self.current_results['reflection'][i, j],
                                    self.current_results['absorption'][i, j]
                                ])
                        
                        writer.writerow(row)
                
                messagebox.showinfo("Export Complete", f"Data exported to {filename}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Error exporting data: {str(e)}")

def main():
    root = tk.Tk()
    app = TMMGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()