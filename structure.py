class StructureManager:
    """Manage different cavity structures"""
    
    def __init__(self):
        pass
    
    def get_structure(self, cavity_type, mirror_material, spacer_material, 
                     active_material, mirror1_thickness, spacer1_thickness,
                     active_thickness, spacer2_thickness, mirror2_thickness):
        """Build structure based on cavity type and materials (legacy method)"""
        
        thicknesses = []
        materials = []
        
        if cavity_type == "full_cavity":
            # Mirror 1 + Spacer 1 + Active + Spacer 2 + Mirror 2
            thicknesses = [
                mirror1_thickness,
                spacer1_thickness,
                active_thickness,
                spacer2_thickness,
                mirror2_thickness
            ]
            materials = [mirror_material, spacer_material, active_material, 
                        spacer_material, mirror_material]
        
        elif cavity_type == "open_cavity":
            # Spacer 1 + Active + Spacer 2
            thicknesses = [
                spacer1_thickness,
                active_thickness,
                spacer2_thickness
            ]
            materials = [spacer_material, active_material, spacer_material]
        
        else:  # no_cavity
            # Just active material
            thicknesses = [active_thickness]
            materials = [active_material]
        
        return thicknesses, materials
    
    def get_flexible_structure_info(self, materials, thicknesses):
        """Get formatted structure information for flexible layer system"""
        info = f"Flexible Layer Structure\n"
        info += f"Number of layers: {len(materials)}\n\n"
        
        for i, (material, thickness) in enumerate(zip(materials, thicknesses)):
            info += f"Layer {i+1}: {material} ({thickness:.1f} nm)\n"
        
        info += f"\nTotal thickness: {sum(thicknesses):.1f} nm\n"
        
        return info
    
    def get_structure_info(self, cavity_type, materials, thicknesses):
        """Get formatted structure information (legacy method)"""
        info = f"Cavity Type: {cavity_type}\n"
        info += f"Materials: {' / '.join(materials)}\n"
        info += f"Thicknesses: {[f'{t:.1f}' for t in thicknesses]} nm\n"
        info += f"Total thickness: {sum(thicknesses):.1f} nm\n"
        
        return info
    
    def validate_structure(self, thicknesses, materials):
        """Validate structure parameters"""
        if len(thicknesses) != len(materials):
            return False, "Number of thicknesses must match number of materials"
        
        if any(t <= 0 for t in thicknesses):
            return False, "All thicknesses must be positive"
        
        if len(materials) == 0:
            return False, "At least one material must be specified"
        
        return True, "Structure is valid"