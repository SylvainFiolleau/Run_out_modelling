from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, QLineEdit, QLabel,
    QTextEdit, QCheckBox, QGroupBox, QToolButton, QFrame, QSizePolicy, QRadioButton, QButtonGroup
)
from PySide6.QtCore import Qt
import geopandas as gpd
import pandas as pd
import sys
import os
import glob
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
print(sys.path)
# --- Dummy imports for illustration, replace with your actual modules ---
import SLBL_Ortho as SLBL
import Prep_Run_Avaframe as prep_Ava
import avaframe
import Consequences as Cons
import wave_start_point as WSP
import functions_Wave as FW


# import WCons

class CollapsibleBox(QWidget):
    def __init__(self, title="", parent=None):
        super().__init__(parent)
        self.toggle_button = QToolButton(text=title, checkable=True, checked=False)
        self.toggle_button.setStyleSheet("QToolButton { border: none; }")
        self.toggle_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle_button.setArrowType(Qt.RightArrow)
        self.toggle_button.clicked.connect(self.on_toggle)
        self.content_area = QFrame()
        self.content_area.setMaximumHeight(0)
        self.content_area.setMinimumHeight(0)
        self.content_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.content_layout = QVBoxLayout()
        self.content_area.setLayout(self.content_layout)
        layout = QVBoxLayout(self)
        layout.addWidget(self.toggle_button)
        layout.addWidget(self.content_area)

    def on_toggle(self):
        checked = self.toggle_button.isChecked()
        self.toggle_button.setArrowType(Qt.DownArrow if checked else Qt.RightArrow)
        self.content_area.setMaximumHeight(16777215 if checked else 0)

    def add_widget(self, widget):
        self.content_layout.addWidget(widget)

    def add_layout(self, layout):
        self.content_layout.addLayout(layout)


class SLBLGui(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SLBL - AvaFrame GUI")
        self.resize(900, 700)
        layout = QVBoxLayout()
        self.path_widgets = {}

        # --- Main Inputs Group ---
        source_box = QGroupBox("Main Inputs")
        source_layout = QVBoxLayout()

        self.add_path_row(source_layout, 'source', "Source Areas File:",
                          "file", "GeoPackage (*.gpkg);;Shapefile (*.shp)", "Select source file")
        self.add_path_row(source_layout, 'dem', "DEM:",
                          "file", "tif (*.tif);;vrt (*.vrt);;All Files (*)", "Select DEM file")

        ids_label_layout = QHBoxLayout()
        self.id_label = QLabel("Select ID(s):")
        self.scenario_label = QLabel("Scenario ID(s):")
        ids_label_layout.addWidget(self.id_label)
        ids_label_layout.addWidget(self.scenario_label)
        source_layout.addLayout(ids_label_layout)

        ids_input_layout = QHBoxLayout()
        self.id_input = QLineEdit()
        self.id_input.setPlaceholderText("e.g., 18001, 18002")
        self.scenario_input = QLineEdit()
        self.scenario_input.setPlaceholderText("If scenario ID selected do not select IDs, e.g., 1800199, 1800101")
        ids_input_layout.addWidget(self.id_input)
        ids_input_layout.addWidget(self.scenario_input)
        source_layout.addLayout(ids_input_layout)

        source_box.setLayout(source_layout)
        layout.addWidget(source_box)

        # --- SLBL Group ---
        SLBL_box = CollapsibleBox("SLBL")
        self.run_slbl = QCheckBox("Run SLBL")
        SLBL_box.add_widget(self.run_slbl)
        self.tolerance_label = QLabel("SLBL Tolerance:")
        self.tol_input = QLineEdit()
        self.tol_input.setPlaceholderText("Tolerance value (e.g., 0.3 or leave blank to select automatically)")
        SLBL_box.add_widget(self.tolerance_label)
        SLBL_box.add_widget(self.tol_input)
        self.add_path_row(SLBL_box.content_layout, 'slbl_out', "SLBL output folder:", "folder")
        layout.addWidget(SLBL_box)

        # --- Avaframe Group ---
        Ava_box = CollapsibleBox("Avaframe Runout Modelling")
        self.run_avaframe = QCheckBox("Run AvaFrame")
        Ava_box.add_widget(self.run_avaframe)
        self.add_path_row(Ava_box.content_layout, 'slbl_out2',
                          "SLBL output folder (same as above):", "folder")
        self.add_path_row(Ava_box.content_layout, 'sim_out',
                          "Simulations output folder:", "folder")
        self.avaframe_autoparam = QCheckBox("Auto parametrization")
        Ava_box.add_widget(self.avaframe_autoparam)
        Ava_box.add_widget(QLabel("Runout Parameters:"))
        self.timeSteps = QLineEdit();
        self.timeSteps.setPlaceholderText("timeSteps (Default: 0.1)")
        self.nb_part = QLineEdit();
        self.nb_part.setPlaceholderText("Number of particles (Default: 10000)")
        self.Resolution = QLineEdit();
        self.Resolution.setPlaceholderText("Resolution (Default: 5)")
        self.fric_label = QLabel("Friction model")
        self.Voellmy = QRadioButton("Voellmy")
        self.Coulomb = QRadioButton("Coulomb")

        self.fric_group = QButtonGroup(self)
        self.fric_group.setExclusive(True)
        self.fric_group.addButton(self.Voellmy)
        self.fric_group.addButton(self.Coulomb)

        fric_hbox = QHBoxLayout()
        fric_hbox.addWidget(self.fric_label)
        fric_hbox.addWidget(self.Voellmy)
        fric_hbox.addWidget(self.Coulomb)

        Ava_box.add_layout(fric_hbox)

        # self.FricModel = QLineEdit(); self.FricModel.setPlaceholderText("Friction Model (Default: Voellmy)")
        self.mu_input = QLineEdit();
        self.mu_input.setPlaceholderText("mu (Default: 0.06)")
        self.xi_input = QLineEdit();
        self.xi_input.setPlaceholderText("turbulence coef. (Default: 400)")
        self.ava_density = QLineEdit();
        self.ava_density.setPlaceholderText("Density (Default: 2700)")
        for widget in [self.timeSteps, self.nb_part, self.Resolution, self.mu_input, self.xi_input, self.ava_density]:
            Ava_box.add_widget(widget)

        # --- Consequences Subgroup ---
        ROCons_box = QGroupBox("Consequences")
        ROCons_layout = QVBoxLayout()
        self.run_consequence = QCheckBox("Runout Consequences")
        ROCons_layout.addWidget(self.run_consequence)
        self.add_path_row(ROCons_layout, 'resident_shp',
                          "Residential shapefile:", "file", "Shapefile (*.shp);;GeoPackage (*.gpkg)",
                          "Select residential shapefile")
        self.add_path_row(ROCons_layout, 'runout_cons_out',
                          "Runout consequences, Output file:", "file",
                          "GeoPackage (*.gpkg);;CSV Files (*.csv);;All Files (*)", "runout_consequences.gpkg",
                          save_dialog=True)
        ROCons_box.setLayout(ROCons_layout)
        Ava_box.add_widget(ROCons_box)
        layout.addWidget(Ava_box)

        # --- Wave Simulation Group ---
        Wave_box = CollapsibleBox("Wave simulation")
        self.run_wave = QCheckBox("Run Wave Simulation")
        Wave_box.add_widget(self.run_wave)
        self.add_path_row(Wave_box.content_layout, 'sim_out2',
                          "Simulations output folder (wave):", "folder")
        self.add_path_row(Wave_box.content_layout, 'start_points',
                          "Wave starting points:", "file", "Shapefile (*.shp);;All Files (*)", "StartPoints.shp",
                          save_dialog=True)
        self.add_path_row(Wave_box.content_layout, 'water_input',
                          "Water area input file:", "file", "tif (*.tif);;All Files (*)", "raster file")
        self.add_path_row(Wave_box.content_layout, 'wave_out',
                          "Wave output folder:", "folder")
        Wave_box.add_widget(QLabel("Wave Parameters:"))
        self.Max_volume = QLineEdit();
        self.Max_volume.setPlaceholderText("Max volume [Mm3] (Default: 50)")
        self.Min_wave = QLineEdit();
        self.Min_wave.setPlaceholderText("Min wave height [m] (Default: 3)")
        self.Max_wave = QLineEdit();
        self.Max_wave.setPlaceholderText("Max wave height [m] (Default: 500)")
        for widget in [self.Max_volume, self.Min_wave, self.Max_wave]:
            Wave_box.add_widget(widget)

        # --- Wave Consequences Subgroup ---
        Cons_box = QGroupBox("Wave consequences")
        Cons_layout = QVBoxLayout()
        self.run_wave_consequence = QCheckBox("Run Wave Consequences")
        Cons_layout.addWidget(self.run_wave_consequence)
        self.add_path_row(Cons_layout, 'resident_shp_wave',
                          "Residential shapefile (wave):", "file", "Shapefile (*.shp);;GeoPackage (*.gpkg)",
                          "Select residential shapefile")
        self.add_path_row(Cons_layout, 'dem_50m',
                          "DEM 50m resolution:", "file", "tif (*.tif);;vrt (*.vrt);;All Files (*)", "Select 50m DEM")
        self.add_path_row(Cons_layout, 'wave_cons_out',
                          "Wave consequences, Output file:", "file",
                          "GeoPackage (*.gpkg);;CSV Files (*.csv);;All Files (*)", "wave_consequences.gpkg",
                          save_dialog=True)
        Cons_box.setLayout(Cons_layout)
        Wave_box.add_widget(Cons_box)
        layout.addWidget(Wave_box)

        # --- Run Button and Log ---
        self.run_button = QPushButton("Run Selected Tasks")
        self.run_button.clicked.connect(self.run_tasks)
        self.output_log = QTextEdit()
        self.output_log.setReadOnly(True)
        layout.addWidget(self.run_button)
        layout.addWidget(self.output_log)
        self.setLayout(layout)

    def add_path_row(self, parent_layout, key, label_text, dialog_type, filter="All Files (*)", placeholder="",
                     save_dialog=False):
        label = QLabel(label_text)
        path_edit = QLineEdit()
        if placeholder:
            path_edit.setPlaceholderText(placeholder)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(lambda: self.browse_path(path_edit, dialog_type, filter, save_dialog))
        row_layout = QHBoxLayout()
        row_layout.addWidget(path_edit)
        row_layout.addWidget(browse_btn)
        parent_layout.addWidget(label)
        parent_layout.addLayout(row_layout)
        self.path_widgets[key] = path_edit

    def browse_path(self, line_edit, dialog_type="file", filter="All Files (*)", save_dialog=False):
        if dialog_type == "file":
            if save_dialog:
                file_path, _ = QFileDialog.getSaveFileName(self, "Select Output File", "", filter)
            else:
                file_path, _ = QFileDialog.getOpenFileName(self, "Select Input File", "", filter)
            if file_path:
                line_edit.setText(file_path)
        elif dialog_type == "folder":
            folder_path = QFileDialog.getExistingDirectory(self, "Select Folder")
            if folder_path:
                line_edit.setText(folder_path)

    def run_tasks(self):
        log = ""
        # --- Gather inputs ---
        source_path_str = self.path_widgets['source'].text()
        dem_path_str = self.path_widgets['dem'].text()
        IDs_selected = self.id_input.text()
        SCENARIOID_selected = self.scenario_input.text()

        IDs = [item.strip() for item in IDs_selected.split(',') if item.strip()]

        SCENARIOID = [item.strip() for item in SCENARIOID_selected.split(',') if item.strip()]
        log += (f"Selected IDs: {IDs}\nSelected ScenarioIDs: {SCENARIOID}\n")

        if not source_path_str:
            log += "No source file selected.\n"
            self.output_log.append(log)
            return
        if not dem_path_str:
            log += "No dem file selected.\n"
            self.output_log.append(log)
            return

        sources_areas = gpd.read_file(source_path_str)
        if IDs:
            print(IDs)
            IDs_int = [int(id_str) for id_str in IDs]
            sources_areas = sources_areas[sources_areas['USTABILEFJELLID'].isin(IDs_int)]
        if SCENARIOID:
            print('Scen')
            SCENARIOID_int = [int(id_str) for id_str in SCENARIOID]

            sources_areas = sources_areas[sources_areas['SCENARIOID'].isin(SCENARIOID_int)]
        print(sources_areas)
        # --- SLBL ---
        if self.run_slbl.isChecked():
            log += "Running SLBL...\n"
            out_folder_SLBL = self.path_widgets['slbl_out'].text().strip()
            resample_resolution = None
            tolerance = self.tol_input.text().strip()
            if not tolerance:
                tolerance = -99
            else:
                try:
                    tolerance = float(tolerance)
                except ValueError:
                    log += "tolerance value not a number\n"
                    self.output_log.append(log)
                    return
            if not out_folder_SLBL:
                log += "No output folder selected.\n"
                self.output_log.append(log)
                return
            SLBL.prepare_and_compute_SLBL(sources_areas, dem_path_str, out_folder_SLBL, tolerance, IDs, SCENARIOID,
                                          resample_resolution)
            log += "Done SLBL.\n"

        # --- AvaFrame ---
        if self.run_avaframe.isChecked():
            ava_folder = r'C:\Users\Fiolleau_Sylvain\Documents\Avaframe_Software\AvaFrame\avaframe'

            if not self.run_slbl.isChecked():
                out_folder_SLBL = self.path_widgets['slbl_out2'].text()
                if not out_folder_SLBL:
                    log += "No SLBL output folder selected.\n"
                    self.output_log.append(log)
                    return
            out_folder = self.path_widgets['sim_out'].text()
            if not out_folder:
                log += "No Simulation output folder selected.\n"
                self.output_log.append(log)
                return
            if self.avaframe_autoparam.isChecked():
                log += f"Runout Params - auto parametrization\n"
                log += "Running AvaFrame ...\n"
                prep_Ava.Prep_Ava_Simul(sources_areas, dem_path_str, ava_folder, out_folder, out_folder_SLBL, IDs,
                                        SCENARIOID, autoparam=True,
                                        ava_time_step=0.1, fric_model='Voellmy', ava_mu_user=0.06, ava_t_coef=400,
                                        nb_part=10000,
                                        out_res=5, ava_density=2700, run_simul=1)
                log += "Done AvaFrame.\n"
            else:
                timeSteps = self.timeSteps.text().strip()
                if not timeSteps:
                    timeSteps = 0.1
                else:
                    try:
                        timeSteps = float(timeSteps)
                    except ValueError:
                        log += "timeSteps value not a number\n"
                        self.output_log.append(log)
                        return
                nb_part = self.nb_part.text().strip()
                if not nb_part:
                    nb_part = 10000
                else:
                    try:
                        nb_part = int(nb_part)
                    except ValueError:
                        log += "nb_part value not an integer\n"
                        self.output_log.append(log)
                        return
                resolution = self.Resolution.text().strip()
                if not resolution:
                    resolution = 5
                else:
                    try:
                        resolution = int(resolution)
                    except ValueError:
                        log += "resolution value not an integer\n"
                        self.output_log.append(log)
                        return
                if self.Voellmy.isChecked():
                    FrictionModel = 'Voellmy'
                elif self.Coulomb.isChecked():
                    FrictionModel = 'Coulomb'
                else:
                    FrictionModel = 'Voellmy'

                mu = self.mu_input.text().strip()
                if not mu:
                    mu = 0.06
                else:
                    try:
                        mu = float(mu)
                    except ValueError:
                        log += "mu value not a number\n"
                        self.output_log.append(log)
                        return
                turbulence = self.xi_input.text().strip()
                if not turbulence:
                    turbulence = 400
                else:
                    try:
                        turbulence = float(turbulence)
                    except ValueError:
                        log += "turbulence value not a number\n"
                        self.output_log.append(log)
                        return
                Density = self.ava_density.text().strip()
                if not Density:
                    Density = 2700
                else:
                    try:
                        Density = float(Density)
                    except ValueError:
                        log += "Density value not a number\n"
                        self.output_log.append(log)
                        return
                log += (
                    f"Runout Params - friction model: {FrictionModel}, mu: {mu}, xi: {turbulence}, density: {Density}\n"
                    f"time steps: {timeSteps}, nb particles: {nb_part}, resolution: {resolution}\n")
                log += "Running AvaFrame ...\n"
                prep_Ava.Prep_Ava_Simul(sources_areas, dem_path_str, ava_folder, out_folder, out_folder_SLBL, IDs,
                                        SCENARIOID, autoparam=False,
                                        ava_time_step=timeSteps, fric_model=FrictionModel, ava_mu_user=mu,
                                        ava_t_coef=turbulence,
                                        nb_part=nb_part, out_res=resolution, ava_density=Density, run_simul=1)
                log += "Done AvaFrame.\n"

        # --- Consequence Analysis ---
        if self.run_consequence.isChecked():
            out_folder = self.path_widgets['sim_out'].text()
            log += "Running Consequence Analysis...\n"
            residential_shp_file = self.path_widgets['resident_shp'].text()
            if not residential_shp_file:
                log += "No residential shapefile selected.\n"
                self.output_log.append(log)
                return
            output_consequence_file_path = self.path_widgets['runout_cons_out'].text()
            if not output_consequence_file_path:
                log += "No output runout consequence file selected.\n"
                self.output_log.append(log)
                return
            Cons.run_out_consequences(out_folder, output_consequence_file_path, residential_shp_file, IDs, SCENARIOID)
            log += "Done Consequence Analysis.\n"

        # --- Wave Simulation ---
        if self.run_wave.isChecked():
            log += "Running Wave Simulation...\n"
            fkb_water_path = self.path_widgets['water_input'].text()
            output_starting_points = self.path_widgets['start_points'].text()
            if not self.run_avaframe.isChecked():
                out_folder = self.path_widgets['sim_out2'].text()
                if not out_folder:
                    log += "No simulation folder selected.\n"
                    self.output_log.append(log)
                    return
            WSP.Launch_wave_start_points(fkb_water_path, out_folder, output_starting_points, IDs, SCENARIOID)
            log += "Done with starting points. Start wave propagation.\n"
            output_dir_wave = self.path_widgets['wave_out'].text()
            if not output_dir_wave:
                log += "No output Wave directory selected.\n"
                self.output_log.append(log)
                return
            max_vol = self.Max_volume.text()
            if not max_vol:
                max_vol = 50
            else:
                try:
                    max_vol = float(max_vol)
                except ValueError:
                    log += "max_vol value not a number\n"
                    self.output_log.append(log)
                    return
            min_wave = self.Min_wave.text()
            if not min_wave:
                min_wave = 3
            else:
                try:
                    min_wave = float(min_wave)
                except ValueError:
                    log += "min_wave value not a number\n"
                    self.output_log.append(log)
                    return
            max_wave = self.Max_wave.text()
            if not max_wave:
                max_wave = 500
            else:
                try:
                    max_wave = float(max_wave)
                except ValueError:
                    log += "max_wave value not a number\n"
                    self.output_log.append(log)
                    return
            lim_azimuth = 100
            parallel_run = 0

            FW.Wave_simulation(output_dir_wave, output_starting_points, max_vol, min_wave, max_wave, False, IDs,
                               SCENARIOID, lim_azimuth, fkb_water_path, parallel_run)

            log += "Done with wave propagation.\n"

        # --- Wave Consequences ---
        if self.run_wave_consequence.isChecked():
            log += "Start wave consequences.\n"
            output_dir_wave = self.path_widgets['wave_out'].text()
            if not output_dir_wave:
                log += "No output Wave directory selected.\n"
                self.output_log.append(log)
                return
            # Example dummy code for demonstration
            # layer_output = gpd.GeoDataFrame(columns=columns_to_keep, geometry="geometry", crs="EPSG:25833")
            dtm50_path = self.path_widgets['dem_50m'].text()
            columns_to_keep = ["USTABILEFJELLID", "SCENARIOID", "Consequences", "geometry"]

            layer_output = gpd.GeoDataFrame(columns=columns_to_keep, geometry="geometry", crs="EPSG:25833")
            wave_consequences_file = self.path_widgets['wave_cons_out'].text()
            if not dtm50_path:
                log += "No output DTM at 50m resolution selected.\n"
                self.output_log.append(log)
                return
            IdsFold = [item for item in IDs if os.path.isdir(os.path.join(output_dir_wave, item))]
            for Id in IdsFold:
                scenarioIds = os.listdir(os.path.join(output_dir_wave, Id))
                for scenarioId in scenarioIds:
                    files = '*.tif'
                    raster_paths = glob.glob(os.path.join(output_dir_wave, Id, scenarioId, files))
                    try:
                        Consequences = WCons.Wave_consequence(raster_paths, residential_shp_file, dtm50_path)
                    except Exception as e:
                        log += f"Error in process_sites {Id}, {scenarioId}: {str(e)}\n"
                        continue
                    if Consequences is not None and not Consequences.empty:
                        layer_output = pd.concat([layer_output, Consequences[columns_to_keep]], ignore_index=True)
            layer_output['geometry'] = layer_output['geometry'].buffer(0)
            layer_output.to_file(wave_consequences_file, driver='GPKG')
            log += "Done wave consequences.\n"

        self.output_log.append(log)


if __name__ == "__main__":
    app = QApplication([])
    gui = SLBLGui()
    gui.show()
    app.exec()
