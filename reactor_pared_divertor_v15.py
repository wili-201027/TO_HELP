# -*- coding: utf-8 -*-
"""
REACTOR_PARED_DIVERTOR_V15.py
============================
Modelo CAD conceptual integrado para Blender 3.6-4.x. Version V15, con exportación de centrolineas para el motor magnetostático.
V14 amplia los canales de refrigeración y rediseña las bobinas superconductoras de integración de auditoría: asientos de soportes, cobertura de capas,
boquilla replegada, modo HV excluido por defecto, validación de parámetros y auditoría
BVH de intersecciones de superficies sólidas.
Incluye primera pared segmentada, blanket conceptual, blindaje, inyector de combustible
con precámara y pistón de resorte, electrodos helicoidales de ensayo HV y acceso de
mantenimiento con cassette/puentes desmontables. La apertura de mantenimiento es un
modo de visualización, no un procedimiento operativo validado.

Incluye:
  - Vasija toroidal con espesor real de 30 mm y pasos de tuberia recortados.
  - Vasija exterior completa por defecto (CUTAWAY_VIEW=False); corte opcional para inspeccionar el interior.
  - Plasma de referencia con los mismos parametros que los divertores.
  - 4 bandas de divertor helicoidal; 16 cassettes por banda (64 en total).
  - Blancos segmentados, respaldos separados en dos alas y ranura de bombeo abierta.
  - Serpentines de refrigeracion en la cara posterior de los blancos.
  - Colectores helicoidales separados de ida/retorno y ramales con rutas en capas.
  - Cuatro conductos del plenum a bombas turbomoleculares, colector de descarga común elevado fuera de la envolvente de bobinas,
    bomba primaria y brida final de conexión a la instalación de vacío. El carril de
    servicio de bombeo se ha retraído para ganar margen frente a la pared de la cámara.
  - Pasos reales en la malla de la vasija para 16 conectores de refrigeracion y 4
    salidas de bombeo. Los conductos y bridas de agua quedan dentro de la envolvente radial
    de 0.50-1.50 m; la estación primaria se eleva sobre el reactor como equipo auxiliar.
  - Soportes de cassettes con anclajes a la cara interior de la vasija.
  - Identificacion visible de conexiones: agua de ida (azul), retorno (rojo) y
    extraccion de vacio (verde), con manguitos/cuerpos de conexion separados.
  - Una ventana de inspeccion/diagnostico con cristal curvado, marco y pernos.
  - Bobinas superconductoras ampliadas con núcleo de cable, aislamiento, carcasa estructural y doble paso interno de helio.
  - Canales separados de suministro/retorno, con puente en U en bobinas cerradas y acoplamientos terminales en segmentos desmontables.
  - Sección exterior aproximada de paquete: 48 mm radial x 64 mm transversal; las tres capas se separan y se permiten salientes fuera del toroide.
  - Dos bastidores toroidales: uno intermedio soporta H/XB y otro exterior soporta XA; puntales cortos sustituyen a los postes que cruzaban bobinas.
  - Sillas, zapatas y postes dimensionados para las nuevas capas de bobinas.
  - Primera pared segmentada, módulos de blanket conceptual y blindaje neutrónico
    en los sectores entre pistas del divertor, con corredores poloidales reservados
    para los cuatro carriles de colectores y huecos alrededor de puertos.
  - Inyector de combustible conceptual: precámara, válvula de llenado, línea de
    evacuación conectada al colector de vacío, pistón accionado por resorte y boquilla.
  - Electrodos de ensayo de alta tensión: dos hélices de sentido opuesto, bornes
    separados y pasamuros/aisladores independientes; sin parámetros eléctricos de operación.
  - Módulos y anclajes de mantenimiento; escotilla de extracción, raíles de servicio,
    cassette de demostración y puentes desmontables en las bobinas que cruzan la escotilla.
  - La envolvente 0.50–1.50 m corresponde solo al toroide/vasija; las bobinas sobresalen deliberadamente.

DIMENSIONES INTERPRETADAS
-------------------------
Los radios 0.50–1.50 m se aplican exclusivamente al toroide/vasija y a su
envolvente toroidal de referencia; no limitan las bobinas ni sus soportes.
Radio mayor del eje toroidal: 1.00 m
Radio menor exterior de vasija: 0.435 m
Espesor nominal de pared: 0.030 m
Radio menor de la cara interior de la vasija: 0.405 m
Las bobinas superconductoras se desplazan a capas externas y pueden sobresalir
de la envolvente 0.50–1.50 m. Su envolvente se informa por separado.
Los canales de agua usan radios de 3.0/3.8/3.2 mm; el paquete de bobina es 48x64 mm
con dos pasajes internos de helio de 2.6 mm de radio.

NOTA FISICA IMPORTANTE
----------------------
Este archivo construye geometria CAD paramétrica. No calcula equilibrio MHD,
separatrices, puntos X, footprints de flujo, carga térmica, tensiones, soldaduras,
fatiga ni hidráulica. La posición helicoidal representa la topología geométrica
solicitada, no demuestra que esas sean las huellas magnéticas reales del plasma.
Los diámetros, claros y materiales son valores de maqueta, no dimensionado de
un reactor operativo.

EJECUCION
---------
Blender > Scripting > Open > reactor_pared_divertor_v14.py > Run Script.
La vasija completa es la vista por defecto. Para inspección en corte, cambiar
CUTAWAY_VIEW=True. Para ver el modo de extracción, cambiar MAINTENANCE_HATCH_OPEN=True.
El objeto legado Plasma_Volume_Combined se conserva por defecto. Solo se elimina si
REMOVE_LEGACY_PLASMA_OBJECT=True; si se conserva, no se crea una segunda malla de plasma.
Para el ensayo HV separado, cambiar HV_TEST_MODE=True; el plasma de referencia se omite.
El modo abierto oculta la tapa, el cassette instalado y los puentes de bobina desmontables;
solo representa una configuración de mantenimiento conceptual.
El script es re-ejecutable: reemplaza su colección anterior sin tocar el resto de la escena. Las rutas funcionales se identifican mediante propiedades propias RDV14_ClearanceRoute/RDV14_ClearanceKind, no por nombre de colección.
"""

import bpy
import bmesh
import math
import json
from mathutils import Vector
from mathutils.kdtree import KDTree
from mathutils.bvhtree import BVHTree

# =============================================================================
# 1. DIMENSIONES PRINCIPALES (metros)
# =============================================================================
R_INNER_AXIS = 0.50       # limite radial minimo de la envolvente del toroide/vasija solamente
R_OUTER_AXIS = 1.50       # limite radial maximo de la envolvente del toroide/vasija solamente
R_MAJOR = 1.000           # radio mayor del eje toroidal
# La vasija se retrae para dejar espacio real a los paquetes de bobinas en el exterior.
R_MINOR_VESSEL_OUTER = 0.435
VESSEL_WALL_THICKNESS = 0.030
R_MINOR_CHAMBER = R_MINOR_VESSEL_OUTER - VESSEL_WALL_THICKNESS  # 0.405 m

# Plasma reducido para reservar volumen a blancos, respaldos, canales, colectores y vacío.
R_PLASMA_BASE = 0.310
AMP_HELICAL = 0.002
AMP_QUADRUPOLAR = 0.0015
PLASMA_MINOR_MAX_BOUND = R_PLASMA_BASE + abs(AMP_HELICAL) + abs(AMP_QUADRUPOLAR)  # .3135 m

# Topología indicada: 4 bobinas helicoidales, 4 vueltas poloidales por vuelta toroidal,
# y 16 bobinas inclinadas (8 pares). El CAD no resuelve el campo de dichas bobinas.
N_HELICAL_COILS = 4
POLOIDAL_TURNS_PER_TOROIDAL_TURN = 4
N_INCLINED_COILS = 16
N_X_PAIRS = N_INCLINED_COILS // 2
N_TRACKS = 4
CASSETTES_PER_TRACK = 16
N_CASSETTES = N_TRACKS * CASSETTES_PER_TRACK
THETA_PHASE = -math.pi / 2.0
TOROIDAL_PITCH = 2.0 * math.pi / CASSETTES_PER_TRACK
CASSETTE_SPAN_FRACTION = 0.88
CASSETTE_SPAN = TOROIDAL_PITCH * CASSETTE_SPAN_FRACTION

# Modo de visualización de la vasija.
CUTAWAY_VIEW = False
CUTAWAY_HALF_ANGLE = math.radians(35.0)  # hueco visual de 70° centrado en phi=0
CREATE_PLASMA_REFERENCE = True
# Por defecto se conserva cualquier objeto legado llamado Plasma_Volume_Combined.
# Activar esta opción solo si se desea eliminarlo explícitamente al regenerar la escena.
REMOVE_LEGACY_PLASMA_OBJECT = False

# Resolución de la vasija con los portales recortados.
VESSEL_N_PHI = 512
VESSEL_N_THETA = 256
CUTAWAY_GUIDE_RADIUS = 0.0012

# =============================================================================
# 2. DIMENSIONADO DE BLANCOS, RESPALDOS, CANALES Y TUBOS (metros/radianes)
# =============================================================================
ROOT_NAME = "Reactor_Pared_Divertor_V15"
PREFIX = "RDV15_"
COLLECTION_NAMES = {
    "plasma": "00_Plasma_Referencia",
    "vessel": "01_Vasija_Con_Pasos",
    "backing": "02_Respaldos_y_Cassettes",
    "targets": "03_Blancos_Tungsteno",
    "baffles": "04_Baffles_y_Pleno",
    "cooling": "05_Canales_Refrigeracion",
    "headers": "06_Colectores_Helicoidales",
    "feedthroughs": "07_Pasos_y_Conectores",
    "supports": "08_Soportes_y_Anclajes",
    "guides": "09_Referencias",
    "coils": "10_Conductores_Bobinas",
    "coil_mounts": "11_Soportes_Bobinas",
    "labels": "12_Etiquetas_Conexiones",
    "window": "13_Ventana_Inspeccion",
    "vacuum_station": "14_Estacion_Bombeo_Vacio",
    "layers": "15_Primera_Pared_Blanket_Blindaje",
    "fuel": "16_Inyeccion_Combustible",
    "ignition": "17_Ignicion_HV_Pruebas",
    "maintenance": "18_Mantenimiento_Cassettes",
    "coil_cooling": "19_Circuitos_Helio_Bobinas",
}

TARGET_HALF_ANGLE = math.radians(24.0)
PUMPING_GAP_ANGLE = math.radians(8.0)
PUMPING_GAP_HALF = PUMPING_GAP_ANGLE / 2.0
BACKING_MARGIN_ANGLE = math.radians(4.0)
BACKING_PHI_MARGIN = 0.0025
TILES_TOROIDAL = 4
TILE_ROWS_PER_WING = 2
TILE_GAP_FRACTION = 0.055

# Superficies anidadas. Los offsets se miden en la dirección poloidal radial de la
# sección toroidal (r_minor), para que placas, tuberías y penetraciones coincidan.
TILE_FACE_OFFSET = 0.014
TILE_THICKNESS = 0.008
TILE_INNER_OFFSET = TILE_FACE_OFFSET
TILE_OUTER_OFFSET = TILE_FACE_OFFSET + TILE_THICKNESS
BACKING_GAP = 0.003
BACKING_INNER_OFFSET = TILE_OUTER_OFFSET + BACKING_GAP  # .025 m
BACKING_THICKNESS = 0.011
BACKING_OUTER_OFFSET = BACKING_INNER_OFFSET + BACKING_THICKNESS  # .036 m
BAFFLE_INNER_OFFSET = 0.010
BAFFLE_OUTER_OFFSET = 0.026
BAFFLE_EDGE_OVERLAP = math.radians(0.25)
BAFFLE_THETA_HALF_WIDTH = math.radians(0.34)
END_FLANGE_PHI = 0.0030

# El centro de los serpentines queda a 2 mm del reverso del respaldo: el tubo
# penetra ~0.3 mm en el respaldo y representa un canal soldado/brazed en contacto.
CHANNEL_CENTER_OFFSET = 0.0380
CHANNEL_RADIUS = 0.0030  # 3.0 mm de radio (+30% respecto a V12)
HEADER_CENTER_OFFSET = 0.0465  # solo referencia de diseño; los colectores V9 usan radios absolutos
HEADER_RADIUS = 0.0038  # 3.8 mm de radio (+27% respecto a V12)
BRANCH_CENTER_OFFSET = 0.0570  # solo referencia; los ramales suben a su carril absoluto
BRANCH_RADIUS = 0.0032  # 3.2 mm de radio (+28% respecto a V12)
# Cuatro carriles de servicio a radios menores absolutos. La separación entre
# carriles evita que ida/retorno dependan de una superficie plasma ondulada.
COOLANT_COLLECTOR_RADII = {
    "A_Ida": 0.354,
    "A_Retorno": 0.365,
    "B_Ida": 0.376,
    "B_Retorno": 0.387,
}
PUMP_SERVICE_MINOR_RADIUS = 0.392
PIPE_CLEARANCE_MIN = 0.001  # mínimo exigido por el comprobador (1 mm)
PIPE_CHECK_MAX_STEP = 0.002  # densificación de chequeo, 2 mm
HEADER_THETA_NEAR = math.radians(34.0)
HEADER_THETA_FAR = math.radians(40.0)
CHANNEL_SAMPLES_PER_PASS = 30
CHANNEL_ROWS = 4

# Un par de soportes por cassette, en los hombros, nunca en la ranura de bombeo.
SUPPORT_RADIUS = 0.0040
SUPPORT_ANCHOR_RADIUS = 0.0050
SUPPORT_WALL_CENTER_MINOR = R_MINOR_CHAMBER - 0.0040
SUPPORT_THETA_OFFSET = TARGET_HALF_ANGLE + math.radians(2.0)

# Penetraciones de refrigeración: un tubo por cada colector helicoidal.
COOLANT_PORT_PHI = math.pi
COOLANT_PORT_HALF_THETA = 0.032
COOLANT_PORT_HALF_PHI = 0.020
COOLANT_STUB_MINOR = 0.492

# Conducto de bombeo por banda. Se coloca en el centro de un cassette, no en una
# junta entre cassettes. El desfase +18° evita que salga directamente por theta=0/pi.
PUMP_PORT_THETA_OFFSET = math.radians(18.0)
PUMP_PORT_HALF_THETA = 0.043
PUMP_PORT_HALF_PHI = 0.028
PUMP_DUCT_OFFSET = 0.0570
PUMP_DUCT_RADIUS = 0.0035
PUMP_DUCT_MOUTH_RADIUS = 0.0060
PUMP_STUB_MINOR = 0.512
# El colector de descargas está sobre el reactor, fuera de la cámara; su ramal conecta
# con una bomba primaria y una brida de salida identificada hacia la instalación de vacío.
# El colector de escape se eleva por fuera de la envolvente de las bobinas X_A.
# El anillo a 0.760 m deja margen a carcasa de bobina (radio conservador 32 mm)
# y a su bastidor exterior (hasta ~0.723 m). La estación primaria se desplaza
# solidariamente para conservar una conexión continua desde el colector.
VACUUM_HEADER_MINOR_RADIUS = 0.760
VACUUM_HEADER_RADIUS = 0.0070
VACUUM_PRIMARY_INLET_MINOR = 0.810
VACUUM_PRIMARY_STACK_SHIFT = VACUUM_PRIMARY_INLET_MINOR - 0.650
VACUUM_EXHAUST_BYPASS_PHI_DEG = {1: 122.5}  # evita X03_A al salir de la turbomolecular T02

# Ventana curvada de inspeccion/diagnostico. El cristal se centra en el espesor
# de pared y el marco queda sobre la cara exterior de la vasija.
WINDOW_PHI = 5.0 * math.pi / 8.0
WINDOW_THETA = math.pi / 4.0
WINDOW_HALF_THETA = math.radians(7.5)
WINDOW_HALF_PHI = math.radians(1.45)
WINDOW_GLASS_MINOR = R_MINOR_VESSEL_OUTER - VESSEL_WALL_THICKNESS * 0.5
WINDOW_FRAME_MINOR = R_MINOR_VESSEL_OUTER + 0.0018
WINDOW_FRAME_RADIUS = 0.0040
WINDOW_BOLT_RADIUS = 0.0030

# Paquetes superconductores: núcleo de cable + aislamiento + carcasa estructural.
# Los radios 0.50–1.50 m SOLO limitan el toroide/vasija. Las bobinas se desplazan
# a capas exteriores; esta envolvente ampliada es intencionada.
COIL_CENTERLINE_MINOR = 0.490       # capa de las 4 bobinas helicoidales
COIL_X_CENTERLINE_MINOR = 0.625     # centro de las dos capas X A/B
COIL_X_SEPARATION_R = 0.040         # A/B separadas 80 mm entre ejes
# Dimensiones exteriores de la carcasa: 48 mm radial x 64 mm transversal.
COIL_CONDUCTOR_RADIUS = 0.024       # semieje radial exterior de carcasa, 24 mm
COIL_TRANSVERSE_HALF_WIDTH = 0.032  # semieje transversal exterior, 32 mm
# Capas internas: núcleo superconductivo, aislamiento y canales de helio en el
# espacio periférico del paquete. Geometría conceptual; no es una sección CICC validada.
COIL_CORE_RADIAL_HALF_WIDTH = 0.0115
COIL_CORE_TRANSVERSE_HALF_WIDTH = 0.0150
COIL_INSULATION_RADIAL_HALF_WIDTH = 0.0230
COIL_INSULATION_TRANSVERSE_HALF_WIDTH = 0.0260
COIL_CASE_RADIAL_HALF_WIDTH = COIL_CONDUCTOR_RADIUS
COIL_CASE_TRANSVERSE_HALF_WIDTH = COIL_TRANSVERSE_HALF_WIDTH
COIL_PROFILE_SEGMENTS = 20
COIL_MIN_SURFACE_CLEARANCE = 0.004  # 4 mm mínimos exigidos entre paquetes exteriores
COIL_COLLISION_CHECK_SAMPLES = 4096
COIL_SADDLE_RADIAL_DEPTH = 0.012
COIL_SUPPORT_TARGET_MINOR = COIL_CENTERLINE_MINOR - COIL_CONDUCTOR_RADIUS - 0.002
COIL_SADDLE_CENTER_MINOR = COIL_CENTERLINE_MINOR - COIL_CONDUCTOR_RADIUS - COIL_SADDLE_RADIAL_DEPTH * 0.5
COIL_SUPPORT_FOOT_MINOR = R_MINOR_VESSEL_OUTER + 0.004
COIL_SUPPORT_FOOT_DEPTH = 0.008
COIL_SUPPORT_POST_START = R_MINOR_VESSEL_OUTER + COIL_SUPPORT_FOOT_DEPTH
COIL_SUPPORT_POST_WIDTH = 0.022
COIL_HELICAL_SADDLE_SIZE = (0.054, 0.064, COIL_SADDLE_RADIAL_DEPTH)
COIL_X_SADDLE_SIZE = (0.054, 0.064, COIL_SADDLE_RADIAL_DEPTH)
COIL_X_TILT = math.radians(20.0)

# V15: referencia de corriente del estudio corregido. Son razones de ensayo,
# no corrientes eléctricas validadas ni una fuente de campo dentro de Blender.
V15_CURRENT_REFERENCE_RATIOS = {
    "H": 1.0,
    "XeqA": -0.955,
    "XeqB": +1.843,
    "XpolA": 0.0,
    "XpolB": 0.0,
}
EXPORT_V15_MAGNETIC_GEOMETRY = True
V15_MAGNETIC_EXPORT_PATH = "//reactor_v15_magnetic_geometry.json"
V15_MAGNETIC_EXPORT_SAMPLES = 4096
COIL_LABEL_MINOR = 0.760
# Bastidores toroidales de soporte intermedio y exterior. Separan físicamente
# cada capa X sin postes que atraviesen bobinas de otras capas.
COIL_INNER_FRAME_MINOR = 0.538
COIL_INNER_FRAME_HALF_THICKNESS = 0.006
COIL_OUTER_FRAME_MINOR = 0.715
COIL_OUTER_FRAME_HALF_THICKNESS = 0.008
COIL_FRAME_LINK_RADIUS = 0.005
# Doble paso interno de helio: canal de suministro y canal de retorno separados,
# unidos por un codo en el extremo opuesto de cada bobina cerrada. Las bobinas
# de mantenimiento abiertas usan un paso continuo con dos interfaces terminales.
COIL_HELIUM_CHANNEL_OFFSET = 0.0200
COIL_HELIUM_CHANNEL_RADIUS = 0.0026
COIL_HELIUM_PORT_STUB_LENGTH = 0.042  # valor nominal; cada familia calcula su alcance al bastidor
COIL_HELIUM_PORT_COLLAR_RADIUS = 0.0065

# Agujeros rectangulares en la placa de fondo del plenum. Sus dimensiones son
# mayores que el conducto para permitir holgura de montaje.
PLENUM_FLOOR_INNER_OFFSET = 0.040
PLENUM_FLOOR_OUTER_OFFSET = 0.043
PLENUM_HOLE_HALF_PHI = 0.022
PLENUM_HOLE_HALF_THETA = 0.032
PLENUM_N_PHI = 512
PLENUM_N_THETA = 20
SLOT_WALL_THICKNESS_ANGLE = math.radians(0.68)
SLOT_WALL_INNER_OFFSET = 0.010
SLOT_WALL_OUTER_OFFSET = 0.041

# --- V14: capas modulares de primera pared/blanket/blindaje (offsets desde plasma) ---
FIRST_WALL_INNER_OFFSET = 0.0390
FIRST_WALL_OUTER_OFFSET = 0.0460
BLANKET_INNER_OFFSET = 0.0490
BLANKET_OUTER_OFFSET = 0.0710
NEUTRON_SHIELD_INNER_OFFSET = 0.0740
NEUTRON_SHIELD_OUTER_OFFSET = 0.0870
GAP_PANEL_THETA_LO = math.radians(29.0)  # comienza fuera del respaldo del divertor (28°)
GAP_PANEL_THETA_HI = math.radians(61.0)  # termina antes del respaldo de la pista siguiente (62°)
# Corredores poloidales abiertos en primera pared/blanket/blindaje para que los
# cuatro colectores helicoidales pasen sin atravesar los paneles. Dos offsets son
# los colectores B de la pista actual; los otros dos corresponden a los colectores
# A de la pista siguiente, expresados en la fase local de esta banda.
COLLECTOR_SERVICE_OFFSETS = tuple(sorted((
    HEADER_THETA_NEAR,
    HEADER_THETA_FAR,
    0.5 * math.pi - HEADER_THETA_FAR,
    0.5 * math.pi - HEADER_THETA_NEAR,
)))
COLLECTOR_SERVICE_HALF_ANGLE = math.radians(0.85)
FW_TILES_TOROIDAL = 4
FW_TILES_POLOIDAL = 4  # referencia de tamaño nominal; la banda se subdivide en intervalos de servicio
BACKING_TILES_TOROIDAL = 4
BACKING_TILES_POLOIDAL = 4
FW_PANEL_PHI_MARGIN = 0.0040
FW_TILE_GAP_FRACTION = 0.055
BACKING_TILE_GAP_FRACTION = 0.055

# Escotilla de servicio alineada con el cassette T01_C02 y un modo desmontable.
MAINTENANCE_TRACK = 0
MAINTENANCE_SEGMENT = 1
MAINTENANCE_HATCH_PHI = (MAINTENANCE_SEGMENT + 0.5) * TOROIDAL_PITCH
MAINTENANCE_HATCH_THETA = (THETA_PHASE + 2.0 * math.pi * MAINTENANCE_TRACK / N_TRACKS
                           + POLOIDAL_TURNS_PER_TOROIDAL_TURN * MAINTENANCE_HATCH_PHI) % (2.0 * math.pi)
MAINTENANCE_HATCH_HALF_THETA = math.radians(29.0)
MAINTENANCE_HATCH_HALF_PHI = 0.5 * CASSETTE_SPAN + math.radians(0.6)
MAINTENANCE_HATCH_OPEN = False  # True = escotilla abierta y cassette de demostración extraído
MAINTENANCE_FRAME_RADIUS = 0.0009
MAINTENANCE_BOLT_RADIUS = 0.0018

# Inyector de prueba: esclusa/precámara exterior, en un claro poloidal entre las pistas.
FUEL_PORT_THETA = math.pi / 2.0
FUEL_PORT_PHI = math.pi / 16.0
FUEL_PORT_HALF_THETA = math.radians(1.0)
FUEL_PORT_HALF_PHI = math.radians(0.55)
FUEL_NOZZLE_TIP_MINOR = 0.3265
FUEL_CHAMBER_START_MINOR = 0.452
FUEL_CHAMBER_END_MINOR = 0.565
FUEL_CHAMBER_OUTER_RADIUS = 0.026
FUEL_CHAMBER_INNER_RADIUS = 0.020
FUEL_PISTON_MINOR = 0.505
FUEL_SPRING_START_MINOR = 0.512
FUEL_SPRING_END_MINOR = 0.552
FUEL_BRANCH_PHI_OFFSET = 0.105
FUEL_NOZZLE_TUBE_RADIUS = 0.0032
FUEL_NOZZLE_MIN_CLEARANCE = 0.0050

# Electrodos de ensayo en modo no operativo: dos hélices opuestas dentro del volumen
# central, retirables/aisladas antes de cualquier régimen de confinamiento.
HV_PORTS = {
    "HV_Anodo": {"theta": math.pi, "phi": math.pi / 3.0, "minor": 0.245, "hand": 1.0},
    "HV_Catodo": {"theta": math.pi, "phi": 5.0 * math.pi / 3.0, "minor": 0.270, "hand": -1.0},
}
HV_TEST_MODE = False  # True: sustituye la escena normal de plasma por la configuracion de ensayo HV
HV_ELECTRODE_PHASE = 5.0 * math.pi / 3.0
HV_ELECTRODE_RADIUS = 0.0018
HV_TERMINAL_MINOR = 0.482

# Colecciones/materiales nombrados para limpiar de forma segura al re-ejecutar.
# Registro de rutas para comprobación integral de tubos contra tubos y bobinas.
CLEARANCE_PATHS = []
RUN_SOLID_MESH_AUDIT = True
SOLID_AUDIT_MAX_OBJECTS = 9000
SOLID_AUDIT_EPSILON = 1e-6
SOLID_AUDIT_MAX_REPORTED = 80
SOLID_AUDIT_MAX_PAIRS = 100000

# =============================================================================
# 3. MATERIALES
# =============================================================================
def material(name, color, metallic=0.4, roughness=0.35, emission=None, alpha=1.0):
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name=name)
    mat.diffuse_color = (color[0], color[1], color[2], alpha)
    mat.metallic = metallic
    mat.roughness = roughness
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    bsdf = nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = (color[0], color[1], color[2], alpha)
        bsdf.inputs["Metallic"].default_value = metallic
        bsdf.inputs["Roughness"].default_value = roughness
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = alpha
        if emission is not None:
            if "Emission Color" in bsdf.inputs:
                bsdf.inputs["Emission Color"].default_value = (*emission[:3], 1.0)
            elif "Emission" in bsdf.inputs:
                bsdf.inputs["Emission"].default_value = (*emission[:3], 1.0)
            if "Emission Strength" in bsdf.inputs:
                bsdf.inputs["Emission Strength"].default_value = emission[3]
    # Compatibilidad de transparencia entre versiones de Blender.
    if alpha < 0.999:
        if hasattr(mat, "surface_render_method"):
            try:
                mat.surface_render_method = "DITHERED"
            except Exception:
                pass
        if hasattr(mat, "blend_method"):
            try:
                mat.blend_method = "HASHED"
            except Exception:
                pass
    return mat


def make_materials():
    return {
        "vessel_inner": material("RDV9_Mat_Vasija_Interior", (0.23, 0.29, 0.34), 0.82, 0.31),
        "vessel_outer": material("RDV9_Mat_Vasija_Exterior", (0.31, 0.37, 0.42), 0.84, 0.29),
        "tungsten": [
            material("RDV9_Mat_Tungsteno_A", (0.43, 0.45, 0.48), 0.86, 0.28),
            material("RDV9_Mat_Tungsteno_B", (0.33, 0.36, 0.40), 0.86, 0.29),
            material("RDV9_Mat_Tungsteno_C", (0.51, 0.53, 0.56), 0.87, 0.25),
            material("RDV9_Mat_Tungsteno_D", (0.39, 0.42, 0.46), 0.84, 0.30),
        ],
        "backing": material("RDV9_Mat_Respaldo", (0.105, 0.135, 0.165), 0.82, 0.31),
        "baffle": material("RDV9_Mat_Baffle", (0.19, 0.235, 0.275), 0.78, 0.34),
        "plenum": material("RDV9_Mat_Pleno_Oscuro", (0.018, 0.027, 0.038), 0.38, 0.42),
        "support": material("RDV9_Mat_Soportes", (0.16, 0.19, 0.22), 0.85, 0.29),
        "channel": material("RDV9_Mat_Canales_Refrigeracion", (0.08, 0.58, 0.72), 0.58, 0.25),
        "supply": material("RDV9_Mat_Refrigeracion_Ida", (0.035, 0.25, 0.83), 0.52, 0.25),
        "return": material("RDV9_Mat_Refrigeracion_Retorno", (0.90, 0.20, 0.055), 0.52, 0.25),
        "pump": material("RDV14_Mat_Conductos_Bombeo", (0.08, 0.55, 0.38), 0.55, 0.28),
        "first_wall": material("RDV14_Mat_Primera_Pared_W", (0.58, 0.60, 0.64), 0.86, 0.24),
        "blanket": material("RDV14_Mat_Blanket_Breeder", (0.58, 0.39, 0.16), 0.78, 0.31),
        "neutron_shield": material("RDV14_Mat_Blindaje_Neutronico", (0.23, 0.34, 0.31), 0.83, 0.34),
        "fuel": material("RDV14_Mat_Linea_Combustible", (0.78, 0.58, 0.12), 0.58, 0.25),
        "fuel_shell": material("RDV14_Mat_Precamara_Inyector", (0.56, 0.40, 0.17), 0.82, 0.27),
        "fuel_piston": material("RDV14_Mat_Piston_Inyector", (0.68, 0.72, 0.76), 0.88, 0.23),
        "fuel_spring": material("RDV14_Mat_Resorte_Inyector", (0.91, 0.64, 0.14), 0.78, 0.23),
        "fuel_gas": material("RDV14_Mat_Gas_Combustible", (0.98, 0.72, 0.10), 0.15, 0.35, emission=(0.8, 0.32, 0.02, 0.5), alpha=0.30),
        "hv_anode": material("RDV14_Mat_Electrodo_Anodo", (0.96, 0.20, 0.12), 0.77, 0.2, emission=(0.8, 0.05, 0.02, 0.25)),
        "hv_cathode": material("RDV14_Mat_Electrodo_Catodo", (0.10, 0.38, 0.98), 0.77, 0.2, emission=(0.02, 0.12, 0.8, 0.25)),
        "hv_insulator": material("RDV14_Mat_Aislante_HV", (0.82, 0.85, 0.90), 0.22, 0.18),
        "maintenance": material("RDV14_Mat_Mantenimiento", (0.76, 0.50, 0.13), 0.80, 0.27),
        "maintenance_latch": material("RDV14_Mat_Cierre_Escotilla", (0.30, 0.34, 0.39), 0.85, 0.25),
        "flange": material("RDV9_Mat_Bridas", (0.20, 0.24, 0.28), 0.85, 0.28),
        "plasma": material("RDV9_Mat_Plasma_Referencia", (0.9, 0.02, 0.35), 0.0, 0.2, emission=(1.0, 0.015, 0.22, 8.0), alpha=0.42),
        "window": material("RDV9_Mat_Cristal_Ventana", (0.10, 0.48, 0.72), 0.12, 0.12, emission=(0.06, 0.28, 0.45, 0.35), alpha=0.28),
        "coil_helix": material("RDV9_Mat_Asiento_Bobina_Helicoidal", (0.035, 0.42, 0.88), 0.78, 0.24),
        "coil_x": material("RDV9_Mat_Asiento_Bobina_Inclinada_X", (0.98, 0.52, 0.08), 0.78, 0.25),
        "coil_helix_conductor": material("RDV9_Mat_Conductor_Helicoidal", (0.035, 0.22, 0.88), 0.91, 0.19),
        "coil_x_eq_a": material("RDV9_Mat_Conductor_X_Ecuatorial_A", (0.035, 0.72, 0.96), 0.91, 0.19),
        "coil_x_eq_b": material("RDV9_Mat_Conductor_X_Ecuatorial_B", (0.10, 0.32, 0.82), 0.91, 0.19),
        "coil_x_pol_a": material("RDV9_Mat_Conductor_X_Poloidal_A", (0.96, 0.24, 0.045), 0.91, 0.19),
        "coil_x_pol_b": material("RDV9_Mat_Conductor_X_Poloidal_B", (1.0, 0.58, 0.06), 0.91, 0.19),
        "coil_case": material("RDV14_Mat_Carcasa_Bobina_Inox", (0.25, 0.31, 0.37), 0.90, 0.24, alpha=0.58),
        "coil_insulation": material("RDV14_Mat_Aislamiento_Bobina", (0.62, 0.78, 0.86), 0.18, 0.22, alpha=0.34),
        "coil_core": material("RDV14_Mat_Nucleo_Superconductor", (0.16, 0.62, 0.88), 0.82, 0.20),
        "coil_helium": material("RDV14_Mat_Canal_Helio", (0.10, 0.88, 0.94), 0.34, 0.18, emission=(0.05, 0.48, 0.62, 0.42)),
        "coil_helium_supply": material("RDV14_Mat_Entrada_Helio", (0.05, 0.38, 1.0), 0.70, 0.18),
        "coil_helium_return": material("RDV14_Mat_Retorno_Helio", (0.70, 0.30, 0.95), 0.70, 0.18),
        "label_water_in": material("RDV9_Mat_Etiqueta_Agua_Ida", (0.035, 0.25, 0.83), 0.35, 0.3),
        "label_water_out": material("RDV9_Mat_Etiqueta_Agua_Retorno", (0.90, 0.20, 0.055), 0.35, 0.3),
        "label_pump": material("RDV9_Mat_Etiqueta_Vacio", (0.06, 0.46, 0.28), 0.35, 0.3),
        "label_text": material("RDV9_Mat_Texto_Claro", (0.92, 0.96, 1.0), 0.1, 0.32, emission=(0.65, 0.75, 0.9, 0.25)),
    }

# =============================================================================
# 4. FUNCIONES DE COORDENADAS
# =============================================================================
def wrap_angle(a):
    return a % (2.0 * math.pi)


def angular_distance(a, b):
    return abs((a - b + math.pi) % (2.0 * math.pi) - math.pi)


def plasma_minor_radius(theta, phi):
    return (
        R_PLASMA_BASE
        + AMP_HELICAL * math.cos(4.0 * theta - 16.0 * phi)
        + AMP_QUADRUPOLAR * math.cos(2.0 * theta - 4.0 * phi)
    )


def torus_point(theta, phi, minor_radius):
    return Vector((
        (R_MAJOR + minor_radius * math.cos(theta)) * math.cos(phi),
        (R_MAJOR + minor_radius * math.cos(theta)) * math.sin(phi),
        minor_radius * math.sin(theta),
    ))


def plasma_point(theta, phi):
    return torus_point(theta, phi, plasma_minor_radius(theta, phi))


def surface_point(theta, phi, offset=0.0):
    """Punto CAD anidado, desplazado por r_minor para alinear los pasos de pared.

    No es el offset diferencial exacto de una superficie helicoidal; es un offset
    paramétrico coherente que conserva theta/phi para que los cassettes, tubos y
    penetraciones coincidan geométricamente.
    """
    return torus_point(theta, phi, plasma_minor_radius(theta, phi) + offset)


def ribbon_theta(phi, theta_offset=0.0, phase=THETA_PHASE):
    return phase + POLOIDAL_TURNS_PER_TOROIDAL_TURN * phi + theta_offset


def phi_for_track_theta(phase, theta_offset, target_theta, sector=0):
    """Escoge una de las cuatro posiciones toroidales equivalentes de una pista helicoidal."""
    base = (target_theta - phase - theta_offset) / POLOIDAL_TURNS_PER_TOROIDAL_TURN
    return wrap_angle(base + (sector % POLOIDAL_TURNS_PER_TOROIDAL_TURN) * (math.pi / 2.0))


def pump_port_phi(track, phase):
    # Las cuatro tomas salen hacia la zona superior (theta=90 grados), en sectores separados.
    return phi_for_track_theta(phase, PUMP_PORT_THETA_OFFSET, math.pi / 2.0, sector=track)


def coolant_header_label(port_name):
    """Convierte T01_A_Ida en A_Ida, para mantener cada salida en su colector."""
    parts = port_name.split("_")
    return "_".join(parts[1:])


def coolant_header_radius(label_or_port_name):
    label = label_or_port_name
    if label not in COOLANT_COLLECTOR_RADII:
        label = coolant_header_label(label_or_port_name)
    return COOLANT_COLLECTOR_RADII[label]


def make_track_point(phi, theta_offset, phase, offset):
    theta = ribbon_theta(phi, theta_offset, phase)
    return surface_point(theta, phi, offset)


def clean_old_system():
    CLEARANCE_PATHS.clear()
    def remove_recursive(col):
        for child in list(col.children):
            remove_recursive(child)
        for obj in list(col.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(col)

    # Sustituir la versión anterior del sistema y evitar que sobrevivan cassettes
    # duplicados de la colección V3 que se entregó antes.
    for name in (ROOT_NAME, "Reactor_Pared_Divertor_V14", "Reactor_Pared_Divertor_V13", "Reactor_Pared_Divertor_V12", "Reactor_Pared_Divertor_V11", "Reactor_Pared_Divertor_V10", "Reactor_Pared_Divertor_V9", "Reactor_Pared_Divertor_V8", "Reactor_Pared_Divertor_V7", "Reactor_Pared_Divertor_V6", "Reactor_Pared_Divertor_V5", "Reactor_Pared_Divertor_V4", "Divertor_3D_Topologico_v3"):
        root = bpy.data.collections.get(name)
        if root is not None:
            remove_recursive(root)
    # El objeto legado no se elimina por defecto: puede pertenecer a otra escena
    # o ser una referencia que el usuario desea conservar. Su borrado es opt-in.
    if REMOVE_LEGACY_PLASMA_OBJECT:
        old_plasma = bpy.data.objects.get("Plasma_Volume_Combined")
        if old_plasma is not None:
            old_mesh = old_plasma.data if old_plasma.type == "MESH" else None
            bpy.data.objects.remove(old_plasma, do_unlink=True)
            if old_mesh is not None and old_mesh.users == 0:
                bpy.data.meshes.remove(old_mesh)
    # Limpiar datos sin usuario que quedaron de ambos sistemas generados.
    for datablocks in (bpy.data.meshes, bpy.data.curves):
        for block in list(datablocks):
            if block.users == 0 and (block.name.startswith((PREFIX, "RDV14_", "RDV13_", "RDV12_", "RDV11_", "RDV10_", "RDV9_", "RDV8_", "RDV7_", "RDV6_", "RDV5_", "RDV4_", "DIV3_V3_"))):
                datablocks.remove(block)


def create_collections():
    root = bpy.data.collections.new(ROOT_NAME)
    bpy.context.scene.collection.children.link(root)
    cols = {"root": root}
    for key, name in COLLECTION_NAMES.items():
        col = bpy.data.collections.new(name)
        root.children.link(col)
        cols[key] = col
    return cols


def link_mesh_object(name, verts, faces, collection, mats=None, smooth=True):
    mesh = bpy.data.meshes.new(PREFIX + name + "_Mesh")
    mesh.from_pydata([tuple(v) for v in verts], [], faces)
    mesh.update()
    # Recalcular normales coherentemente para superficies cerradas con puertos.
    try:
        bm = bmesh.new()
        bm.from_mesh(mesh)
        if bm.faces:
            bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])
        bm.to_mesh(mesh)
        bm.free()
    except Exception:
        pass
    if smooth:
        for poly in mesh.polygons:
            poly.use_smooth = True
    obj = bpy.data.objects.new(PREFIX + name, mesh)
    collection.objects.link(obj)
    if mats:
        if not isinstance(mats, (list, tuple)):
            mats = [mats]
        for mat in mats:
            mesh.materials.append(mat)
    return obj


def create_curve(name, points, collection, mat, radius, cyclic=False, point_radii=None, resolution=1, route_kind=None, assembly_id=None):
    if len(points) < 2:
        return None
    curve = bpy.data.curves.new(PREFIX + name + "_Curve", type="CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = resolution
    curve.bevel_depth = radius
    curve.bevel_resolution = 3
    curve.fill_mode = "FULL"
    spline = curve.splines.new("POLY")
    spline.points.add(len(points) - 1)
    for i, p in enumerate(points):
        spline.points[i].co = (p.x, p.y, p.z, 1.0)
        if point_radii is not None:
            spline.points[i].radius = point_radii[i]
    spline.use_cyclic_u = cyclic
    obj = bpy.data.objects.new(PREFIX + name, curve)
    collection.objects.link(obj)
    if mat:
        obj.data.materials.append(mat)
    # El registro funcional depende de una marca explícita del objeto, nunca del
    # nombre de la colección. Así se pueden mover objetos entre colecciones sin
    # perder su clasificación en el comprobador de holguras.
    if route_kind is not None:
        valid_route_kinds = {"coil", "cooling", "coolant", "vacuum", "fuel", "hv", "cryo_internal"}
        if route_kind not in valid_route_kinds:
            raise ValueError(f"route_kind no reconocido para {name}: {route_kind}")
        scales = list(point_radii) if point_radii is not None else [1.0] * len(points)
        if len(scales) != len(points):
            raise ValueError(f"La escala variable de {name} no coincide con el número de puntos.")
        # La sección elíptica de las bobinas se envuelve conservadoramente con
        # una esfera de radio igual a la semianchura mayor para evitar falsos libres.
        collision_radius = max(radius, COIL_TRANSVERSE_HALF_WIDTH) if route_kind == "coil" else radius
        obj["RDV14_ClearanceRoute"] = True
        obj["RDV14_ClearanceKind"] = route_kind
        if assembly_id is not None:
            obj["RDV14_CoilAssembly"] = str(assembly_id)
        obj["RDV14_ClearanceRadius"] = float(collision_radius)
        obj["RDV14_ClearanceBaseRadius"] = float(radius)
        obj["RDV14_ClearancePointCount"] = len(points)
        CLEARANCE_PATHS.append({
            "name": name, "object_name": obj.name, "kind": route_kind, "radius": collision_radius,
            "assembly_id": str(assembly_id) if assembly_id is not None else None,
            "scales": scales, "points": [Vector(p).copy() for p in points],
        })
    return obj


def create_uv_sphere(name, center, radius, collection, mat, segments=12, rings=8):
    verts = [Vector(center) + Vector((0, 0, radius))]
    for ri in range(1, rings):
        theta = math.pi * ri / rings
        for si in range(segments):
            phi = 2 * math.pi * si / segments
            verts.append(Vector(center) + Vector((radius * math.sin(theta) * math.cos(phi), radius * math.sin(theta) * math.sin(phi), radius * math.cos(theta))))
    bottom = len(verts)
    verts.append(Vector(center) + Vector((0, 0, -radius)))
    faces = []
    for j in range(segments):
        faces.append((0, 1 + j, 1 + (j + 1) % segments))
    for ring in range(rings - 2):
        a0 = 1 + ring * segments
        b0 = a0 + segments
        for j in range(segments):
            faces.append((a0 + j, b0 + j, b0 + (j + 1) % segments, a0 + (j + 1) % segments))
    last = 1 + (rings - 2) * segments
    for j in range(segments):
        faces.append((last + j, bottom, last + (j + 1) % segments))
    obj = link_mesh_object(name, verts, faces, collection, mat)
    obj.location = (0, 0, 0)
    return obj

# =============================================================================
# 5. MALLAS SÓLIDAS DE PARCHE Y CINTAS HELICOIDALES
# =============================================================================
def create_solid_patch(name, phi0, phi1, theta_off0, theta_off1, phase,
                       inner_offset, outer_offset, collection, mat,
                       n_phi=8, n_theta=4, tile_gaps=None):
    """Parche sólido cerrado entre dos offsets de superficie."""
    verts = []
    # layer 0 = cara hacia plasma; layer 1 = cara posterior.
    for off in (inner_offset, outer_offset):
        for i in range(n_phi + 1):
            phi = phi0 + (phi1 - phi0) * i / n_phi
            for j in range(n_theta + 1):
                theta_off = theta_off0 + (theta_off1 - theta_off0) * j / n_theta
                verts.append(make_track_point(phi, theta_off, phase, off))
    layer_size = (n_phi + 1) * (n_theta + 1)
    row = n_theta + 1
    def ix(layer, i, j):
        return layer * layer_size + i * row + j
    faces = []
    for i in range(n_phi):
        for j in range(n_theta):
            # Separación opcional de pequeñas teselas mediante una malla de celdas.
            if tile_gaps is not None and tile_gaps(i, j, n_phi, n_theta):
                continue
            # Cara posterior orientada hacia la pared.
            faces.append((ix(1, i, j), ix(1, i + 1, j), ix(1, i + 1, j + 1), ix(1, i, j + 1)))
            # Cara interior orientada hacia el plasma.
            faces.append((ix(0, i, j + 1), ix(0, i + 1, j + 1), ix(0, i + 1, j), ix(0, i, j)))
    # Perímetros de cierre entre las dos caras.
    for i in range(n_phi):
        faces.append((ix(0, i, 0), ix(1, i, 0), ix(1, i + 1, 0), ix(0, i + 1, 0)))
        faces.append((ix(0, i + 1, n_theta), ix(1, i + 1, n_theta), ix(1, i, n_theta), ix(0, i, n_theta)))
    for j in range(n_theta):
        faces.append((ix(0, 0, j + 1), ix(1, 0, j + 1), ix(1, 0, j), ix(0, 0, j)))
        faces.append((ix(0, n_phi, j), ix(1, n_phi, j), ix(1, n_phi, j + 1), ix(0, n_phi, j + 1)))
    return link_mesh_object(name, verts, faces, collection, mat)


def create_cyclic_solid_ribbon(name, phase, theta_off0, theta_off1, inner_offset, outer_offset,
                               collection, mat, n_phi=PLENUM_N_PHI, n_theta=PLENUM_N_THETA,
                               hole=None, holes=None):
    """Cinta helicoidal cerrada y sólida. Puede llevar una abertura rectangular."""
    verts = []
    nphi = n_phi
    ntheta = n_theta
    layer_size = nphi * (ntheta + 1)
    row = ntheta + 1
    for off in (inner_offset, outer_offset):
        for i in range(nphi):
            phi = 2.0 * math.pi * i / nphi
            for j in range(ntheta + 1):
                theta_off = theta_off0 + (theta_off1 - theta_off0) * j / ntheta
                verts.append(make_track_point(phi, theta_off, phase, off))
    def ix(layer, i, j):
        return layer * layer_size + (i % nphi) * row + j
    # Celdas omitidas en los agujeros de extracción o mantenimiento.
    removed = [[False for _ in range(ntheta)] for _ in range(nphi)]
    all_holes = []
    if hole is not None:
        all_holes.append(hole)
    if holes:
        all_holes.extend(list(holes))
    for hphi, htheta_off, halfphi, halftheta in all_holes:
        for i in range(nphi):
            pc = 2.0 * math.pi * (i + 0.5) / nphi
            dp = angular_distance(pc, hphi)
            if dp > halfphi:
                continue
            for j in range(ntheta):
                tc = theta_off0 + (theta_off1 - theta_off0) * (j + 0.5) / ntheta
                if abs(tc - htheta_off) <= halftheta:
                    removed[i][j] = True
    faces = []
    for i in range(nphi):
        ni = (i + 1) % nphi
        for j in range(ntheta):
            if removed[i][j]:
                continue
            # Cara exterior a la pared.
            faces.append((ix(1, i, j), ix(1, ni, j), ix(1, ni, j + 1), ix(1, i, j + 1)))
            # Cara hacia el plasma.
            faces.append((ix(0, i, j + 1), ix(0, ni, j + 1), ix(0, ni, j), ix(0, i, j)))
    # Cierres longitudinales de los dos bordes poloidales.
    for i in range(nphi):
        ni = (i + 1) % nphi
        faces.append((ix(0, i, 0), ix(1, i, 0), ix(1, ni, 0), ix(0, ni, 0)))
        faces.append((ix(0, ni, ntheta), ix(1, ni, ntheta), ix(1, i, ntheta), ix(0, i, ntheta)))
    # Cerrar paredes laterales del agujero y cualquier borde de celda omitida.
    for i in range(nphi):
        ni = (i + 1) % nphi
        for j in range(ntheta):
            if not removed[i][j]:
                continue
            # Vecino phi-: borde en phi=i
            im = (i - 1) % nphi
            if not removed[im][j]:
                faces.append((ix(1, i, j), ix(1, i, j + 1), ix(0, i, j + 1), ix(0, i, j)))
            # Vecino phi+: borde en phi=i+1
            ip = (i + 1) % nphi
            if not removed[ip][j]:
                faces.append((ix(1, ip, j + 1), ix(1, ip, j), ix(0, ip, j), ix(0, ip, j + 1)))
            # Vecino theta-: borde j
            if j == 0 or not removed[i][j - 1]:
                faces.append((ix(1, i, j), ix(1, ni, j), ix(0, ni, j), ix(0, i, j)))
            # Vecino theta+: borde j+1
            if j == ntheta - 1 or not removed[i][j + 1]:
                faces.append((ix(1, ni, j + 1), ix(1, i, j + 1), ix(0, i, j + 1), ix(0, ni, j + 1)))
    return link_mesh_object(name, verts, faces, collection, mat)

# =============================================================================
# 6. VASIJA TOROIDAL CON AGUJEROS REALES DE PENETRACIÓN
# =============================================================================
def build_port_specs():
    ports = []
    # Tomas individuales para ida/retorno de ambas alas. Se escogen sectores que
    # dejan las penetraciones lejos del eje central y fuera de la zona de plasma.
    definitions = (
        (-HEADER_THETA_FAR, "A_Ida", 0.0),
        (-HEADER_THETA_NEAR, "A_Retorno", 0.0),
        ( HEADER_THETA_NEAR, "B_Ida", 0.0),
        ( HEADER_THETA_FAR, "B_Retorno", 0.0),
    )
    for track in range(N_TRACKS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_TRACKS
        for idx, (off, label, target_theta) in enumerate(definitions):
            phi = phi_for_track_theta(phase, off, target_theta, sector=(track + idx) % 4)
            theta = wrap_angle(ribbon_theta(phi, off, phase))
            ports.append({
                "name": f"T{track+1:02d}_{label}", "kind": "coolant",
                "theta": theta, "phi": phi,
                "half_theta": COOLANT_PORT_HALF_THETA,
                "half_phi": COOLANT_PORT_HALF_PHI,
                "minor_stub": COOLANT_STUB_MINOR,
            })
    # Un conducto de extraccion por cada ranura helicoidal; la salida se reubica
    # en una de cuatro posiciones equivalentes con el eje aproximadamente vertical.
    for track in range(N_TRACKS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_TRACKS
        phi = pump_port_phi(track, phase)
        theta = wrap_angle(ribbon_theta(phi, PUMP_PORT_THETA_OFFSET, phase))
        ports.append({
            "name": f"T{track+1:02d}_Bombeo", "kind": "pump",
            "theta": theta, "phi": phi,
            "half_theta": PUMP_PORT_HALF_THETA,
            "half_phi": PUMP_PORT_HALF_PHI,
            "minor_stub": PUMP_STUB_MINOR,
        })
    # Precámara de inyección en un claro poloidal, con pasamuros propio.
    ports.append({
        "name": "Precámara_Inyector_Combustible", "kind": "fuel",
        "theta": FUEL_PORT_THETA, "phi": FUEL_PORT_PHI,
        "half_theta": FUEL_PORT_HALF_THETA,
        "half_phi": FUEL_PORT_HALF_PHI,
        "minor_stub": 0.485,
    })
    # Bornes HV opcionales: se omiten por completo en la configuración normal.
    if HV_TEST_MODE:
        for name, spec in HV_PORTS.items():
            ports.append({
                "name": name, "kind": "hv",
                "theta": spec["theta"], "phi": spec["phi"],
                "half_theta": math.radians(1.10),
                "half_phi": math.radians(0.62),
                "minor_stub": HV_TERMINAL_MINOR,
            })
    ports.append({
        "name": "Escotilla_Mantenimiento_Cassette_T01_C02", "kind": "maintenance",
        "theta": MAINTENANCE_HATCH_THETA, "phi": MAINTENANCE_HATCH_PHI,
        "half_theta": MAINTENANCE_HATCH_HALF_THETA,
        "half_phi": MAINTENANCE_HATCH_HALF_PHI,
        "minor_stub": R_MINOR_VESSEL_OUTER,
    })
    ports.append({
        "name": "Ventana_Inspeccion", "kind": "window",
        "theta": WINDOW_THETA, "phi": WINDOW_PHI,
        "half_theta": WINDOW_HALF_THETA,
        "half_phi": WINDOW_HALF_PHI,
        "minor_stub": R_MINOR_VESSEL_OUTER,
    })
    return ports

def port_cell_removed(phi_center, theta_center, ports):
    for port in ports:
        if (angular_distance(phi_center, port["phi"]) <= port["half_phi"]
                and angular_distance(theta_center, port["theta"]) <= port["half_theta"]):
            return True
    return False


def create_vessel_shell(collection, mats, ports):
    """Malla cerrada de la vasija: caras interior/exterior y túneles en los portales."""
    nphi, ntheta = VESSEL_N_PHI, VESSEL_N_THETA
    grid_size = nphi * ntheta
    verts = []
    # Capa exterior, luego capa interior.
    for minor in (R_MINOR_VESSEL_OUTER, R_MINOR_CHAMBER):
        for i in range(nphi):
            phi = 2.0 * math.pi * i / nphi
            for j in range(ntheta):
                theta = 2.0 * math.pi * j / ntheta
                verts.append(torus_point(theta, phi, minor))
    def ix(layer, i, j):
        return layer * grid_size + (i % nphi) * ntheta + (j % ntheta)
    removed = [[False for _ in range(ntheta)] for _ in range(nphi)]
    for i in range(nphi):
        phi_c = 2.0 * math.pi * (i + 0.5) / nphi
        for j in range(ntheta):
            theta_c = 2.0 * math.pi * (j + 0.5) / ntheta
            cutaway = CUTAWAY_VIEW and angular_distance(phi_c, 0.0) < CUTAWAY_HALF_ANGLE
            removed[i][j] = cutaway or port_cell_removed(phi_c, theta_c, ports)
    faces = []
    # Dos superficies anidadas. Las caras eliminadas definen aperturas reales.
    for i in range(nphi):
        ni = (i + 1) % nphi
        for j in range(ntheta):
            nj = (j + 1) % ntheta
            if removed[i][j]:
                continue
            # Exterior: normal que sale de la vasija.
            faces.append((ix(0, i, j), ix(0, ni, j), ix(0, ni, nj), ix(0, i, nj)))
            # Interior: cara invertida, orientada hacia la cámara.
            faces.append((ix(1, i, nj), ix(1, ni, nj), ix(1, ni, j), ix(1, i, j)))
    # Cierre del espesor en cada contorno de apertura (portales y corte visual).
    # Los cantos llevan el material exterior para que se vea el espesor de pared.
    for i in range(nphi):
        ni = (i + 1) % nphi
        for j in range(ntheta):
            if not removed[i][j]:
                continue
            nj = (j + 1) % ntheta
            im = (i - 1) % nphi
            ip = (i + 1) % nphi
            jm = (j - 1) % ntheta
            jp = (j + 1) % ntheta
            if not removed[im][j]:
                faces.append((ix(0, i, j), ix(0, i, nj), ix(1, i, nj), ix(1, i, j)))
            if not removed[ip][j]:
                faces.append((ix(0, ip, nj), ix(0, ip, j), ix(1, ip, j), ix(1, ip, nj)))
            if not removed[i][jm]:
                faces.append((ix(0, i, j), ix(0, ni, j), ix(1, ni, j), ix(1, i, j)))
            if not removed[i][jp]:
                faces.append((ix(0, ni, jp), ix(0, i, jp), ix(1, i, jp), ix(1, ni, jp)))
    obj = link_mesh_object("Vasija_Toroidal_30mm_con_Puertos", verts, faces, collection,
                           [mats["vessel_outer"], mats["vessel_inner"]])
    # Clasificar por radio menor del centro de cada cara. Es más robusto que
    # depender del orden de polígonos tras la reconstrucción de normales con bmesh.
    for poly in obj.data.polygons:
        c = poly.center
        rho = math.sqrt(c.x * c.x + c.y * c.y)
        minor = math.sqrt((rho - R_MAJOR) ** 2 + c.z * c.z)
        poly.material_index = 1 if abs(minor - R_MINOR_CHAMBER) < abs(minor - R_MINOR_VESSEL_OUTER) else 0
    obj["R_major_m"] = R_MAJOR
    obj["Radio_menor_exterior_m"] = R_MINOR_VESSEL_OUTER
    obj["Radio_menor_interior_m"] = R_MINOR_CHAMBER
    obj["Espesor_nominal_m"] = VESSEL_WALL_THICKNESS
    obj["Puertos_pasantes"] = sum(p["kind"] != "window" for p in ports)
    obj["Ventanas"] = sum(p["kind"] == "window" for p in ports)
    obj["Modo_corte"] = bool(CUTAWAY_VIEW)
    return obj


def create_port_collar(name, port, collection, mats):
    theta0, phi0 = port["theta"], port["phi"]
    ht, hp = port["half_theta"] * 1.12, port["half_phi"] * 1.12
    pts = []
    for i in range(96):
        a = 2.0 * math.pi * i / 96
        theta = theta0 + ht * math.cos(a)
        phi = phi0 + hp * math.sin(a)
        pts.append(torus_point(theta, phi, R_MINOR_VESSEL_OUTER + 0.0012))
    return create_curve(name, pts, collection, mats["flange"], 0.0018 if port["kind"] == "coolant" else 0.0024, cyclic=True)


def create_port_tube(port, collection, mats):
    """Manguito recto que pasa por el espesor de la vasija hasta el colector exterior."""
    theta, phi = port["theta"], port["phi"]
    inner_minor = R_MINOR_CHAMBER
    outer_minor = port["minor_stub"]
    start_minor = coolant_header_radius(port["name"]) if port["kind"] == "coolant" else PUMP_SERVICE_MINOR_RADIUS
    # El primer punto coincide exactamente con la línea central del colector; después
    # el manguito atraviesa el espacio de cámara, el espesor y termina fuera.
    radii = []
    points = []
    steps = 20
    for i in range(steps + 1):
        f = i / steps
        r = start_minor + (outer_minor - start_minor) * f
        points.append(torus_point(theta, phi, r))
        radii.append(1.0)
    if port["kind"] == "coolant":
        rad, mat = HEADER_RADIUS, mats["supply"] if "Ida" in port["name"] else mats["return"]
    else:
        rad, mat = PUMP_DUCT_RADIUS, mats["pump"]
    return create_curve("Paso_" + port["name"], points, collection, mat, rad, cyclic=False, point_radii=radii,
                         route_kind="coolant" if port["kind"] == "coolant" else "vacuum")

# =============================================================================
# 6B. VENTANA, CONEXIONES IDENTIFICADAS Y SOPORTES DE BOBINAS
# =============================================================================
def torus_basis(theta, phi):
    """Base tangencial local (toroidal, poloidal, normal) de la vasija."""
    tangent_phi = Vector((-math.sin(phi), math.cos(phi), 0.0)).normalized()
    tangent_theta = Vector((-math.sin(theta) * math.cos(phi),
                            -math.sin(theta) * math.sin(phi),
                            math.cos(theta))).normalized()
    normal = Vector((math.cos(theta) * math.cos(phi),
                     math.cos(theta) * math.sin(phi),
                     math.sin(theta))).normalized()
    return tangent_phi, tangent_theta, normal


def create_oriented_box(name, theta, phi, minor_center, dimensions, collection, mat):
    """Bloque orientado localmente a la superficie toroidal; dimensions=(phi,theta,normal)."""
    center = torus_point(theta, phi, minor_center)
    e1, e2, e3 = torus_basis(theta, phi)
    hx, hy, hz = (d * 0.5 for d in dimensions)
    verts = []
    for z in (-hz, hz):
        for y in (-hy, hy):
            for x in (-hx, hx):
                verts.append(center + e1*x + e2*y + e3*z)
    # Indices: lower layer 0..3; upper layer 4..7
    faces = [
        (0, 1, 3, 2), (4, 6, 7, 5),
        (0, 4, 5, 1), (1, 5, 7, 3),
        (3, 7, 6, 2), (2, 6, 4, 0),
    ]
    return link_mesh_object(name, verts, faces, collection, mat, smooth=False)


def create_oriented_cylinder(name, theta, phi, minor_center, radius, depth, collection, mat, segments=24):
    """Cilindro con el eje normal a la superficie exterior de la vasija."""
    center = torus_point(theta, phi, minor_center)
    e1, e2, normal = torus_basis(theta, phi)
    verts = []
    for h in (-depth*0.5, depth*0.5):
        c = center + normal*h
        for i in range(segments):
            a = 2.0*math.pi*i/segments
            verts.append(c + radius*(math.cos(a)*e1 + math.sin(a)*e2))
    faces = []
    for i in range(segments):
        j=(i+1)%segments
        faces.append((i, j, segments+j, segments+i))
    faces.append(tuple(range(segments-1, -1, -1)))
    faces.append(tuple(range(segments, 2*segments)))
    return link_mesh_object(name, verts, faces, collection, mat, smooth=True)


def create_toroidal_support_shell(name, minor_center, half_thickness, collection, mat,
                                 frame_id, n_phi=192, n_theta=64):
    """Bastidor toroidal cerrado como dos superficies concéntricas de pared gruesa.

    El intervalo [minor_center-half_thickness, minor_center+half_thickness]
    representa la sección estructural del aro. No es un cálculo de tensiones.
    """
    radii = (minor_center - half_thickness, minor_center + half_thickness)
    verts, faces = [], []
    ring_size = n_phi * n_theta
    for minor in radii:
        for i in range(n_phi):
            phi = 2.0 * math.pi * i / n_phi
            for j in range(n_theta):
                theta = 2.0 * math.pi * j / n_theta
                verts.append(torus_point(theta, phi, minor))
    for side in range(2):
        base = side * ring_size
        for i in range(n_phi):
            ni = (i + 1) % n_phi
            for j in range(n_theta):
                nj = (j + 1) % n_theta
                a = base + i*n_theta + j
                b = base + ni*n_theta + j
                c = base + ni*n_theta + nj
                d = base + i*n_theta + nj
                faces.append((a, b, c, d) if side == 1 else (d, c, b, a))
    obj = link_mesh_object(name, verts, faces, collection, mat, smooth=True)
    obj["RDV14_SupportFrameId"] = str(frame_id)
    obj["RDV14_SupportFrameShell"] = True
    obj["Radio_menor_centro_m"] = float(minor_center)
    obj["Semiespesor_m"] = float(half_thickness)
    obj["Funcion"] = "Bastidor estructural toroidal conceptual de soporte de bobinas"
    return obj


def create_coil_frame_link(name, theta, phi, start_minor, end_minor,
                           frame_id, assembly_id, cols, mats):
    """Puntal corto normal a la superficie entre una bobina y su bastidor."""
    depth = abs(end_minor - start_minor)
    if depth <= 1e-5:
        return None
    center_minor = 0.5 * (start_minor + end_minor)
    obj = create_oriented_cylinder(
        name, theta, phi, center_minor, COIL_FRAME_LINK_RADIUS,
        depth + 0.002, cols["coil_mounts"], mats["support"], segments=20,
    )
    obj["RDV14_SupportFrameId"] = str(frame_id)
    obj["RDV14_CoilAssembly"] = str(assembly_id)
    obj["Funcion"] = "Puntal de carga entre paquete superconductivo y bastidor toroidal"
    return obj


def create_surface_text(name, text, theta, phi, minor, collection, mat, size=0.014):
    """Texto orientado hacia fuera, sobre una pequeña placa o una superficie curva."""
    font = bpy.data.curves.new(PREFIX + name + "_Font", type="FONT")
    font.body = text
    font.size = size
    font.align_x = "CENTER"
    font.align_y = "CENTER"
    font.extrude = 0.00015
    font.bevel_depth = 0.00005
    obj = bpy.data.objects.new(PREFIX + name, font)
    collection.objects.link(obj)
    obj.location = torus_point(theta, phi, minor)
    obj.rotation_mode = "QUATERNION"
    obj.rotation_quaternion = torus_basis(theta, phi)[2].to_track_quat("Z", "Y")
    if mat:
        obj.data.materials.append(mat)
    return obj


def create_window_assembly(port, cols, mats):
    """Ventana curvada y sellada dentro de la abertura recortada de la vasija."""
    theta0, phi0 = port["theta"], port["phi"]
    ht, hp = port["half_theta"]*0.91, port["half_phi"]*0.91
    ntheta, nphi = 24, 16
    verts, faces = [], []
    for i in range(nphi+1):
        phi = phi0 - hp + 2.0*hp*i/nphi
        for j in range(ntheta+1):
            theta = theta0 - ht + 2.0*ht*j/ntheta
            verts.append(torus_point(theta, phi, WINDOW_GLASS_MINOR))
    stride = ntheta+1
    for i in range(nphi):
        for j in range(ntheta):
            a=i*stride+j; b=(i+1)*stride+j
            faces.append((a, b, b+1, a+1))
    glass = link_mesh_object("Cristal_Ventana_Inspeccion", verts, faces, cols["window"], mats["window"], smooth=True)
    glass["Tipo"] = "Cristal diagnostico conceptual; requiere obturador y validacion radiologica"
    # Marco eliptico curvo, ligeramente por fuera de la vasija, y pernos de sujecion.
    frame_pts=[]
    for i in range(97):
        a=2.0*math.pi*i/96
        th=theta0 + port["half_theta"]*1.08*math.cos(a)
        ph=phi0 + port["half_phi"]*1.08*math.sin(a)
        frame_pts.append(torus_point(th, ph, WINDOW_FRAME_MINOR))
    create_curve("Marco_Metalico_Ventana", frame_pts, cols["window"], mats["flange"], WINDOW_FRAME_RADIUS, cyclic=True, resolution=1)
    for i in range(8):
        a=2.0*math.pi*i/8
        th=theta0 + port["half_theta"]*1.20*math.cos(a)
        ph=phi0 + port["half_phi"]*1.20*math.sin(a)
        create_uv_sphere(f"Perno_Ventana_{i+1:02d}", torus_point(th, ph, WINDOW_FRAME_MINOR+0.001), WINDOW_BOLT_RADIUS, cols["window"], mats["flange"], segments=12, rings=8)
    # Etiqueta desplazada para no tapar el propio cristal.
    create_surface_text("Etiqueta_Ventana", "VENTANA / DIAGNOSTICO", theta0, phi0+0.09, 0.525,
                        cols["labels"], mats["label_text"], size=0.014)


def create_annular_port_flange(name, theta, phi, minor_center, outer_radius, inner_radius,
                               depth, collection, mat, segments=32):
    """Brida anular con taladro central para no taponar el conducto de agua."""
    center = torus_point(theta, phi, minor_center)
    e1, e2, normal = torus_basis(theta, phi)
    verts = []
    # Anillos: exterior/inferior, interior/inferior, exterior/superior, interior/superior.
    for h, radius in ((-depth*0.5, outer_radius), (-depth*0.5, inner_radius),
                      ( depth*0.5, outer_radius), ( depth*0.5, inner_radius)):
        c = center + normal*h
        for i in range(segments):
            a = 2.0*math.pi*i/segments
            verts.append(c + radius*(math.cos(a)*e1 + math.sin(a)*e2))
    def ix(layer, i):
        return layer*segments + (i % segments)
    faces = []
    for i in range(segments):
        j = (i+1) % segments
        # Caras planas de ambos lados de la arandela.
        faces.append((ix(2,i), ix(2,j), ix(3,j), ix(3,i)))
        faces.append((ix(0,j), ix(0,i), ix(1,i), ix(1,j)))
        # Superficie exterior e interior del espesor.
        faces.append((ix(0,i), ix(0,j), ix(2,j), ix(2,i)))
        faces.append((ix(1,j), ix(1,i), ix(3,i), ix(3,j)))
    obj = link_mesh_object(name, verts, faces, collection, mat, smooth=True)
    obj["Tipo"] = "Brida anular con paso central libre"
    return obj


def create_coolant_connector_and_label(port, cols, mats):
    """Brida exterior compacta e identificador agua ida/retorno dentro de la envolvente."""
    theta, phi = port["theta"], port["phi"]
    is_supply = "Ida" in port["name"]
    mat = mats["supply"] if is_supply else mats["return"]
    label_mat = mats["label_water_in"] if is_supply else mats["label_water_out"]
    stub = port["minor_stub"]
    # La brida termina antes del final del tubo; deja libre su sección de paso.
    flange_minor = stub - 0.004
    create_annular_port_flange("Brida_Agua_"+port["name"], theta, phi, flange_minor,
                               outer_radius=0.010, inner_radius=0.0034, depth=0.008,
                               collection=cols["feedthroughs"], mat=mat, segments=32)
    # Placa desplazada en dirección poloidal para que no tape la boca del tubo.
    # La anchura de 60 mm permite identificar conexiones incluso en pares cercanos.
    # El desplazamiento +4.5° crea separación vertical respecto al eje de la toma (theta=0).
    tag_minor = 0.496
    tag_theta = wrap_angle(theta + math.radians(4.5))
    create_oriented_box("Placa_ID_"+port["name"], tag_theta, phi, tag_minor,
                        (0.060, 0.026, 0.004), cols["labels"], label_mat)
    label = ("AGUA IDA " if is_supply else "AGUA RET. ") + port["name"].split("_")[0]
    create_surface_text("Texto_ID_"+port["name"], label, tag_theta, phi, tag_minor+0.0015,
                        cols["labels"], mats["label_text"], size=0.0090)


def create_vacuum_pump_module(track, phase, cols, mats):
    """Bomba turbomolecular conceptual: entrada conectada al conducto del plenum y descarga al colector comun."""
    phi = pump_port_phi(track, phase)
    theta = wrap_angle(ribbon_theta(phi, PUMP_PORT_THETA_OFFSET, phase))
    # Conector entre el pasamuros exterior y la boca de aspiracion de la bomba.
    create_curve(f"T{track+1:02d}_Enlace_Bomba_Vacio",
                 [torus_point(theta, phi, PUMP_STUB_MINOR), torus_point(theta, phi, 0.535)],
                 cols["feedthroughs"], mats["pump"], 0.0065, resolution=1, route_kind="vacuum")
    # Cuerpo de bomba: el eje coincide con la normal local y su boca inferior toca el enlace.
    create_oriented_cylinder(f"T{track+1:02d}_Carcasa_Bomba_Turbomolecular", theta, phi, 0.545,
                             0.018, 0.022, cols["vacuum_station"], mats["pump"], segments=32)
    create_oriented_cylinder(f"T{track+1:02d}_Motor_Bomba_Turbomolecular", theta, phi, 0.564,
                             0.012, 0.016, cols["vacuum_station"], mats["flange"], segments=24)
    # Descarga al colector elevado. T02 necesita un bypass toroidal a r=0.572 m:
    # su trayectoria radial original coincidía con el paquete X03_A. Primero se
    # desplaza por el pasillo libre hasta 122.5°, después asciende fuera de las
    # envolventes de las bobinas. Las otras tres salidas conservan su azimut.
    target_phi = wrap_angle(math.radians(VACUUM_EXHAUST_BYPASS_PHI_DEG.get(track, math.degrees(phi))))
    delta_phi = (target_phi - phi + math.pi) % (2.0 * math.pi) - math.pi
    escape_points = []
    arc_steps = max(1, int(math.ceil(abs(delta_phi) / math.radians(0.20))))
    for k in range(arc_steps + 1):
        f = k / arc_steps
        escape_points.append(torus_point(theta, phi + delta_phi * f, 0.572))
    for k in range(1, 161):
        f = k / 160.0
        minor = 0.572 + (VACUUM_HEADER_MINOR_RADIUS - 0.572) * f
        escape_points.append(torus_point(theta, target_phi, minor))
    escape_obj = create_curve(f"T{track+1:02d}_Escape_Bomba",
                 escape_points, cols["feedthroughs"], mats["pump"], 0.0060,
                 resolution=1, route_kind="vacuum")
    if escape_obj:
        escape_obj["Bypass_azimutal_grados"] = math.degrees(delta_phi)
        escape_obj["Destino"] = "Colector anular de vacío elevado a r_minor=0.760 m"
    # Etiqueta desplazada 7° alrededor del toro para no bloquear el tubo de descarga vertical.
    tag_minor = 0.585
    tag_phi = wrap_angle(phi + math.radians(7.0 if track % 2 == 0 else -7.0))
    create_oriented_box(f"T{track+1:02d}_Placa_ID_Vacio", theta, tag_phi, tag_minor,
                        (0.070, 0.024, 0.004), cols["labels"], mats["label_pump"])
    create_surface_text(f"T{track+1:02d}_Texto_ID_Vacio", f"TURBO T{track+1:02d}", theta, tag_phi,
                        tag_minor+0.003, cols["labels"], mats["label_text"], size=0.0095)


def build_vacuum_exhaust_system(cols, mats):
    """Conecta las cuatro descargas a un colector anular y a una bomba primaria explícita."""
    theta = math.pi / 2.0
    # Colector horizontal sobre el reactor: los cuatro tubos de descarga desembocan aquí.
    samples = 512
    ring_points = [torus_point(theta, 2.0 * math.pi * i / samples, VACUUM_HEADER_MINOR_RADIUS)
                   for i in range(samples)]
    create_curve("Colector_Exhaust_Vacio", ring_points, cols["feedthroughs"],
                 mats["pump"], VACUUM_HEADER_RADIUS, cyclic=True, resolution=1, route_kind="vacuum")

    # Derivación desde el colector hasta la boca de aspiración de la bomba primaria.
    phi_primary = 0.0
    create_curve("Linea_Principal_Vacio",
                 [torus_point(theta, phi_primary, VACUUM_HEADER_MINOR_RADIUS),
                  torus_point(theta, phi_primary, VACUUM_PRIMARY_INLET_MINOR)],
                 cols["feedthroughs"], mats["pump"], 0.0090, resolution=1, route_kind="vacuum")

    # Bancada y cuerpo de la bomba primaria desplazados para que su boca inferior
    # coincida con el nuevo colector a 0.760 m (inlet a 0.810 m).
    shift = VACUUM_PRIMARY_STACK_SHIFT
    create_oriented_box("Bancada_Bomba_Primaria_Vacio", theta, phi_primary, 0.654 + shift,
                        (0.150, 0.115, 0.008), cols["vacuum_station"], mats["support"])
    create_oriented_cylinder("Carcasa_Bomba_Primaria_Vacio", theta, phi_primary, 0.670 + shift,
                             0.030, 0.040, cols["vacuum_station"], mats["pump"], segments=40)
    # Acoplamiento superior y motor eléctrico.
    create_oriented_cylinder("Acoplamiento_Bomba_Primaria", theta, phi_primary, 0.696 + shift,
                             0.016, 0.012, cols["vacuum_station"], mats["flange"], segments=28)
    create_oriented_cylinder("Motor_Electrico_Bomba_Primaria", theta, phi_primary, 0.7245 + shift,
                             0.023, 0.045, cols["vacuum_station"], mats["support"], segments=32)
    # La descarga no queda abierta en el aire: termina en una brida de interfaz a la línea de planta.
    create_curve("Linea_Escape_Primaria",
                 [torus_point(theta, phi_primary, 0.747 + shift), torus_point(theta, phi_primary, 0.780 + shift)],
                 cols["feedthroughs"], mats["pump"], 0.0070, resolution=1, route_kind="vacuum")
    create_oriented_cylinder("Brida_Salida_A_Sistema_Vacio", theta, phi_primary, 0.785 + shift,
                             0.014, 0.010, cols["vacuum_station"], mats["flange"], segments=32)
    create_oriented_box("Placa_ID_Bomba_Primaria", theta, phi_primary, 0.805 + shift,
                        (0.170, 0.050, 0.004), cols["labels"], mats["label_pump"])
    create_surface_text("Texto_ID_Bomba_Primaria", "BOMBA PRIMARIA / ESCAPE A PLANTA",
                        theta, phi_primary, 0.808 + shift, cols["labels"], mats["label_text"], size=0.010)


def create_coil_support_post(name, theta, phi, coil_center_minor, coil_radius,
                             saddle_size, coil_mat, cols, mats):
    """Zapata reforzada + poste + silla para el paquete superconductivo ampliado.

    Todos los apoyos usan una zapata de 8 mm y un poste de sección 22x22 mm.
    Las bobinas X exteriores requieren postes más largos; la comprobación de
    holguras evalúa esos postes contra los demás paquetes de bobina.
    """
    is_helical = name.startswith("H")
    foot_depth = COIL_SUPPORT_FOOT_DEPTH
    foot_minor = R_MINOR_VESSEL_OUTER + foot_depth * 0.5
    post_start = R_MINOR_VESSEL_OUTER + foot_depth

    create_oriented_box(name + "_Zapata", theta, phi, foot_minor,
                        (0.042, 0.050, foot_depth), cols["coil_mounts"], mats["support"])

    contact_minor = coil_center_minor - coil_radius
    saddle_center = contact_minor - saddle_size[2] * 0.5
    saddle_inner = contact_minor - saddle_size[2]
    post_depth = saddle_inner - post_start
    if post_depth > 1e-5:
        post_center = post_start + post_depth * 0.5
        create_oriented_box(name + "_Poste", theta, phi, post_center,
                            (COIL_SUPPORT_POST_WIDTH, COIL_SUPPORT_POST_WIDTH, post_depth),
                            cols["coil_mounts"], mats["support"])
    elif post_depth < -1e-5:
        raise ValueError(f"Soporte {name}: el asiento invade la zapata; profundidad de poste={post_depth:.6f} m")

    create_oriented_box(name + "_Silla", theta, phi, saddle_center, saddle_size,
                        cols["coil_mounts"], coil_mat)


def coil_local_frame(points, index, cyclic=True):
    """Marco local radial/transversal; usa tangente unilateral en bobinas abiertas."""
    n = len(points)
    p = points[index]
    if cyclic:
        tangent = (points[(index + 1) % n] - points[(index - 1) % n]).normalized()
    elif index == 0:
        tangent = (points[1] - points[0]).normalized()
    elif index == n - 1:
        tangent = (points[-1] - points[-2]).normalized()
    else:
        tangent = (points[index + 1] - points[index - 1]).normalized()
    rho = math.sqrt(p.x*p.x + p.y*p.y)
    minor = math.sqrt((rho - R_MAJOR)**2 + p.z*p.z)
    if minor < 1e-9 or rho < 1e-9:
        normal = Vector((1.0, 0.0, 0.0))
    else:
        ctheta = (rho - R_MAJOR) / minor
        stheta = p.z / minor
        normal = Vector((ctheta * p.x / rho, ctheta * p.y / rho, stheta))
    # Elimina pequeñas desviaciones numéricas de la ortogonalidad con la tangente.
    normal = (normal - tangent * normal.dot(tangent)).normalized()
    binormal = tangent.cross(normal).normalized()
    return normal, binormal


def create_elliptical_coil_pack(name, points, collection, mat,
                                radial_half_width=COIL_CONDUCTOR_RADIUS,
                                transverse_half_width=COIL_TRANSVERSE_HALF_WIDTH,
                                profile_segments=COIL_PROFILE_SEGMENTS,
                                cyclic=True, cap_ends=True):
    """Crea la envolvente sólida elíptica 48x64 mm siguiendo el eje de bobina.

    El semieje radial es el menor para respetar el apilamiento de las tres capas;
    el ancho mayor va en la dirección transversal tangente a la superficie toroidal.
    """
    n = len(points)
    verts = []
    ring = profile_segments
    frames = [coil_local_frame(points, i, cyclic=cyclic) for i in range(n)]
    for i, p in enumerate(points):
        normal, binormal = frames[i]
        for j in range(ring):
            a = 2.0 * math.pi * j / ring
            v = p + normal * (radial_half_width * math.cos(a)) + binormal * (transverse_half_width * math.sin(a))
            verts.append(v)
    faces = []
    span_count = n if cyclic else n - 1
    for i in range(span_count):
        ni = (i + 1) % n
        for j in range(ring):
            nj = (j + 1) % ring
            faces.append((i*ring+j, ni*ring+j, ni*ring+nj, i*ring+nj))
    if not cyclic and cap_ends and n >= 2:
        faces.append(tuple(0 * ring + j for j in range(ring)))
        faces.append(tuple((n - 1) * ring + j for j in reversed(range(ring))))
    obj = link_mesh_object(name, verts, faces, collection, mat, smooth=True)
    obj["Seccion_radial_mm"] = 2.0 * radial_half_width * 1000.0
    obj["Ancho_transversal_mm"] = 2.0 * transverse_half_width * 1000.0
    obj["Geometria"] = "Paquete eliptico barrido; eje menor radial para conservar las holguras"
    return obj


def create_superconducting_coil_assembly(name, points, collection, coil_mat, cols, mats,
                                         assembly_id, cyclic=True):
    """Paquete superconductivo multicapa con doble paso interno de helio.

    En bobinas cerradas: dos semipasos paralelos (suministro/retorno) unidos por
    un codo en U en el extremo opuesto; ambos acoplamientos externos quedan juntos
    en el terminal de servicio. En segmentos desmontables abiertos se dibuja un
    paso continuo y dos interfaces terminales para permitir puentes flexibles.

    La geometría es conceptual: no representa un cable CICC validado ni define
    caudal, presión, temperatura, estabilidad criogénica o impregnación eléctrica.
    """
    objects = []
    case = create_elliptical_coil_pack(
        name + "_Carcasa_Estructural_48x64mm", points, collection, mats["coil_case"],
        radial_half_width=COIL_CASE_RADIAL_HALF_WIDTH,
        transverse_half_width=COIL_CASE_TRANSVERSE_HALF_WIDTH,
        cyclic=cyclic,
    )
    insulation = create_elliptical_coil_pack(
        name + "_Aislamiento_Electrico_y_Canales", points, collection, mats["coil_insulation"],
        radial_half_width=COIL_INSULATION_RADIAL_HALF_WIDTH,
        transverse_half_width=COIL_INSULATION_TRANSVERSE_HALF_WIDTH,
        cyclic=cyclic,
    )
    core = create_elliptical_coil_pack(
        name + "_Haz_Superconductor", points, collection, coil_mat,
        radial_half_width=COIL_CORE_RADIAL_HALF_WIDTH,
        transverse_half_width=COIL_CORE_TRANSVERSE_HALF_WIDTH,
        cyclic=cyclic,
    )
    for obj, layer in ((case, "Carcasa_presion_estructural"),
                       (insulation, "Aislamiento_con_pasos_de_helio"),
                       (core, "Haz_superconductor")):
        obj["RDV14_CoilAssembly"] = str(assembly_id)
        obj["RDV14_CoilLayer"] = layer
        obj["Refrigeracion"] = "Doble paso de helio; representacion CAD conceptual"
        objects.append(obj)

    if len(points) >= 8:
        frames = [coil_local_frame(points, i, cyclic=cyclic) for i in range(len(points))]

        def offset_point(index, sign):
            return points[index] + frames[index][1] * (sign * COIL_HELIUM_CHANNEL_OFFSET)

        if cyclic:
            # Suministro: primera mitad del recorrido.
            midpoint = len(points) // 2
            supply_indices = list(range(0, midpoint + 1))
            # Retorno: segunda mitad, recorrida en sentido inverso hasta el mismo terminal.
            return_indices = [0] + list(range(len(points) - 1, midpoint - 1, -1))
            supply_points = [offset_point(i, +1.0) for i in supply_indices]
            return_points = [offset_point(i, -1.0) for i in return_indices]

            supply = create_curve(
                name + "_Canal_Helio_Suministro", supply_points, cols["coil_cooling"],
                mats["coil_helium_supply"], COIL_HELIUM_CHANNEL_RADIUS,
                cyclic=False, resolution=1, route_kind="cryo_internal", assembly_id=assembly_id,
            )
            ret = create_curve(
                name + "_Canal_Helio_Retorno", return_points, cols["coil_cooling"],
                mats["coil_helium_return"], COIL_HELIUM_CHANNEL_RADIUS,
                cyclic=False, resolution=1, route_kind="cryo_internal", assembly_id=assembly_id,
            )
            # Codo de retorno en U alrededor del extremo opuesto del núcleo.
            mid_index = midpoint
            arc_points = []
            for k in range(17):
                a = math.pi * k / 16.0
                arc_points.append(
                    points[mid_index]
                    + frames[mid_index][1] * (COIL_HELIUM_CHANNEL_OFFSET * math.cos(a))
                    + frames[mid_index][0] * (COIL_HELIUM_CHANNEL_OFFSET * math.sin(a))
                )
            jumper = create_curve(
                name + "_Puente_U_Helio_Extremo", arc_points, cols["coil_cooling"],
                mats["coil_helium"], COIL_HELIUM_CHANNEL_RADIUS,
                cyclic=False, resolution=1, route_kind="cryo_internal", assembly_id=assembly_id,
            )
            if supply:
                supply["Circuito"] = "Helio_entrada"
            if ret:
                ret["Circuito"] = "Helio_retorno"
            if jumper:
                jumper["Funcion"] = "Union interna suministro-retorno en extremo opuesto"

            port_specs = (
                (offset_point(0, +1.0), frames[0][0], "Suministro", mats["coil_helium_supply"]),
                (offset_point(0, -1.0), frames[0][0], "Retorno", mats["coil_helium_return"]),
            )
        else:
            # Segmento desmontable: un paso continuo con interfaces a ambos extremos.
            channel_points = [offset_point(i, +1.0) for i in range(len(points))]
            channel = create_curve(
                name + "_Canal_Helio_Segmento", channel_points, cols["coil_cooling"],
                mats["coil_helium"], COIL_HELIUM_CHANNEL_RADIUS,
                cyclic=False, resolution=1, route_kind="cryo_internal", assembly_id=assembly_id,
            )
            if channel:
                channel["Circuito"] = "Segmento desmontable; conectar mediante latiguillos criogenicos"
            port_specs = (
                (offset_point(0, +1.0), frames[0][0], "Suministro_Extremo_A", mats["coil_helium_supply"]),
                (offset_point(len(points)-1, +1.0), frames[-1][0], "Retorno_Extremo_B", mats["coil_helium_return"]),
            )

        # Cada familia lleva los extremos del circuito hasta el bastidor toroidal
        # correspondiente: H y X_B al aro intermedio; X_A al aro exterior.
        if str(assembly_id).startswith("H"):
            frame_id = "INTERMEDIO"
            port_sign = +1.0
            target_minor = COIL_INNER_FRAME_MINOR - COIL_INNER_FRAME_HALF_THICKNESS
        elif str(assembly_id).endswith("_B"):
            frame_id = "INTERMEDIO"
            port_sign = -1.0
            target_minor = COIL_INNER_FRAME_MINOR + COIL_INNER_FRAME_HALF_THICKNESS
        else:
            frame_id = "EXTERIOR"
            port_sign = +1.0
            target_minor = COIL_OUTER_FRAME_MINOR - COIL_OUTER_FRAME_HALF_THICKNESS

        # Longitud calculada según el radio menor real del eje de esta bobina.
        if str(assembly_id).startswith("H"):
            center_minor = COIL_CENTERLINE_MINOR
        elif str(assembly_id).endswith("_B"):
            center_minor = COIL_X_CENTERLINE_MINOR - COIL_X_SEPARATION_R
        else:
            center_minor = COIL_X_CENTERLINE_MINOR + COIL_X_SEPARATION_R
        port_length = abs(target_minor - center_minor)

        for port_point, normal, label, mat in port_specs:
            endpoint = port_point + normal * (port_sign * port_length)
            stub = create_curve(
                name + f"_Toma_Helio_{label}",
                [port_point, port_point.lerp(endpoint, 0.35), endpoint],
                cols["coil_cooling"], mat, COIL_HELIUM_CHANNEL_RADIUS,
                cyclic=False, resolution=1, route_kind="cryo_internal", assembly_id=assembly_id,
            )
            if stub:
                stub["Refrigerante"] = "Helio"
                stub["Tipo_toma"] = label
                stub["Conectada_a_canal_interno"] = True
                stub["RDV14_SupportFrameId"] = frame_id
            collar = create_uv_sphere(
                name + f"_Acoplamiento_Helio_{label}", endpoint,
                COIL_HELIUM_PORT_COLLAR_RADIUS, cols["coil_cooling"], mat,
                segments=20, rings=12,
            )
            collar["RDV14_CoilAssembly"] = str(assembly_id)
            collar["RDV14_SupportFrameId"] = frame_id
            collar["Interfaz_criogenica"] = label
            collar["Destino"] = f"Interfaz al bastidor de distribución de helio {frame_id}; crioplanta no incluida"
            objects.append(collar)
    return objects


def ellipse_support_radius(direction, normal, binormal):
    """Radio soporte de la seccion eliptica en una direccion dada."""
    u = direction.normalized() if direction.length > 1e-12 else Vector((1.0, 0.0, 0.0))
    dn = u.dot(normal)
    db = u.dot(binormal)
    return math.sqrt((COIL_CONDUCTOR_RADIUS * dn)**2 + (COIL_TRANSVERSE_HALF_WIDTH * db)**2)


def build_coil_centerline_samples(sample_count=COIL_COLLISION_CHECK_SAMPLES):
    """Genera muestras de las 20 líneas centrales en la geometría final usada por Blender."""
    specs = []
    for track in range(N_HELICAL_COILS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_HELICAL_COILS
        points = []
        for i in range(sample_count):
            phi = 2.0 * math.pi * i / sample_count
            theta = phase + POLOIDAL_TURNS_PER_TOROIDAL_TURN * phi
            points.append(torus_point(theta, phi, COIL_CENTERLINE_MINOR))
        specs.append((f"H{track+1:02d}", points, COIL_CONDUCTOR_RADIUS))

    for pos in range(8):
        phi_center = pos * (2.0 * math.pi / 8.0)
        is_equatorial = (pos % 2 == 0)
        for coil_index, direction in (("A", 1.0), ("B", -1.0)):
            r_actual = COIL_X_CENTERLINE_MINOR + (
                COIL_X_SEPARATION_R if coil_index == "A" else -COIL_X_SEPARATION_R
            )
            points = []
            for j in range(sample_count):
                theta = 2.0 * math.pi * j / sample_count
                if is_equatorial:
                    phi_offset = direction * COIL_X_TILT * math.sin(theta)
                else:
                    phi_offset = direction * COIL_X_TILT * math.cos(theta)
                phi = phi_center + phi_offset
                points.append(torus_point(theta, phi, r_actual))
            specs.append((f"X{pos+1:02d}{coil_index}", points, COIL_CONDUCTOR_RADIUS))
    return specs


def check_coil_clearances():
    """Chequea las secciones elipticas, no solo las lineas centrales circulares."""
    specs = build_coil_centerline_samples()
    trees, max_steps, frames = [], [], []
    for _name, points, _radius in specs:
        tree = KDTree(len(points))
        max_step = 0.0
        for i, point in enumerate(points):
            tree.insert(point, i)
            d = (points[(i + 1) % len(points)] - point).length
            max_step = max(max_step, d)
        tree.balance()
        trees.append(tree)
        max_steps.append(max_step)
        frames.append([coil_local_frame(points, k) for k in range(len(points))])

    worst_gap = float("inf")
    worst_pair = ("", "")
    clashes = []
    for i, (name_a, points_a, _radius_a) in enumerate(specs):
        for j in range(i + 1, len(specs)):
            name_b, points_b, _radius_b = specs[j]
            min_gap = float("inf")
            # Barrido A -> B. La cota descuenta el error de muestreo de ambas curvas.
            for ia, point_a in enumerate(points_a):
                co, ib, distance = trees[j].find(point_a)
                if distance <= 1e-12:
                    gap = -COIL_CONDUCTOR_RADIUS * 2.0
                else:
                    direction = point_a - co
                    na, ba = frames[i][ia]
                    nb, bb = frames[j][ib]
                    support_a = ellipse_support_radius(direction, na, ba)
                    support_b = ellipse_support_radius(-direction, nb, bb)
                    gap = distance - support_a - support_b - 0.5*(max_steps[i] + max_steps[j])
                if gap < min_gap:
                    min_gap = gap
            # Barrido inverso para no depender de la fase discreta del muestreo.
            for ib, point_b in enumerate(points_b):
                co, ia, distance = trees[i].find(point_b)
                if distance <= 1e-12:
                    gap = -COIL_CONDUCTOR_RADIUS * 2.0
                else:
                    direction = point_b - co
                    nb, bb = frames[j][ib]
                    na, ba = frames[i][ia]
                    support_b = ellipse_support_radius(direction, nb, bb)
                    support_a = ellipse_support_radius(-direction, na, ba)
                    gap = distance - support_a - support_b - 0.5*(max_steps[i] + max_steps[j])
                if gap < min_gap:
                    min_gap = gap
            if min_gap < worst_gap:
                worst_gap, worst_pair = min_gap, (name_a, name_b)
            if min_gap < COIL_MIN_SURFACE_CLEARANCE:
                clashes.append((name_a, name_b, min_gap))

    # Se comprueba la envolvente de los postes frente al perfil eliptico de las otras bobinas.
    support_half_diagonal = math.sqrt(2.0) * COIL_SUPPORT_POST_WIDTH * 0.5
    worst_support_gap = float("inf")
    worst_support = ("", "")
    support_clashes = []
    support_sample_count = 17

    def check_post(support_name, target_name, stheta, sphi, saddle_center, saddle_size, is_helical):
        nonlocal worst_support_gap, worst_support
        foot_depth = COIL_SUPPORT_FOOT_DEPTH
        post_start = R_MINOR_VESSEL_OUTER + foot_depth
        saddle_bottom = saddle_center - saddle_size[2] * 0.5
        post_length = max(0.0, saddle_bottom - post_start)
        post_step = post_length / max(1, support_sample_count - 1)
        for sample in range(support_sample_count):
            r_minor = post_start + post_length * sample / (support_sample_count - 1)
            point = torus_point(stheta, sphi, r_minor)
            for j, (coil_name, _coil_points, _coil_radius) in enumerate(specs):
                if coil_name == target_name:
                    continue
                co, idx, distance = trees[j].find(point)
                direction = point - co
                nframe, bframe = frames[j][idx]
                coil_support = ellipse_support_radius(direction, nframe, bframe)
                gap = distance - 0.5*(max_steps[j] + post_step) - support_half_diagonal - coil_support
                if gap < worst_support_gap:
                    worst_support_gap, worst_support = gap, (support_name, coil_name)
                if gap < COIL_MIN_SURFACE_CLEARANCE:
                    support_clashes.append((support_name, coil_name, gap))

    for track in range(N_HELICAL_COILS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_HELICAL_COILS
        for seg in range(CASSETTES_PER_TRACK):
            phi = 2.0 * math.pi * (seg + 0.5) / CASSETTES_PER_TRACK
            theta = wrap_angle(phase + POLOIDAL_TURNS_PER_TOROIDAL_TURN * phi)
            saddle_center = COIL_CENTERLINE_MINOR - COIL_CONDUCTOR_RADIUS - COIL_HELICAL_SADDLE_SIZE[2] * 0.5
            if MAINTENANCE_HATCH_OPEN and track == MAINTENANCE_TRACK and seg == MAINTENANCE_SEGMENT:
                continue
            check_post(f"H{track+1:02d}_S{seg+1:02d}", f"H{track+1:02d}", theta, phi,
                       saddle_center, COIL_HELICAL_SADDLE_SIZE, True)

    # Los soportes X ya no son postes desde la vasija: son enlaces cortos entre
    # las bobinas y dos bastidores anulares. Se validan sus holguras radiales globales.
    h_outer = math.sqrt((COIL_CENTERLINE_MINOR + COIL_CONDUCTOR_RADIUS)**2 + COIL_TRANSVERSE_HALF_WIDTH**2)
    xb_center = COIL_X_CENTERLINE_MINOR - COIL_X_SEPARATION_R
    xb_inner = math.sqrt((xb_center - COIL_CONDUCTOR_RADIUS)**2 + COIL_TRANSVERSE_HALF_WIDTH**2)
    xa_center = COIL_X_CENTERLINE_MINOR + COIL_X_SEPARATION_R
    xa_outer = math.sqrt((xa_center + COIL_CONDUCTOR_RADIUS)**2 + COIL_TRANSVERSE_HALF_WIDTH**2)
    frame_gaps = {
        "H_a_bastidor_intermedio": (COIL_INNER_FRAME_MINOR - COIL_INNER_FRAME_HALF_THICKNESS) - h_outer,
        "B_a_bastidor_intermedio": xb_inner - (COIL_INNER_FRAME_MINOR + COIL_INNER_FRAME_HALF_THICKNESS),
        "A_a_bastidor_exterior": (COIL_OUTER_FRAME_MINOR - COIL_OUTER_FRAME_HALF_THICKNESS) - xa_outer,
    }
    for label, gap in frame_gaps.items():
        if gap < COIL_MIN_SURFACE_CLEARANCE:
            support_clashes.append((label, "bastidor/paquete", gap))
    worst_frame_gap = min(frame_gaps.values())

    print("\n========== CONTROL DE INTERFERENCIAS DE BOBINAS V14 ==========")
    print(f"Conductores comprobados: {len(specs)} (4 helicoidales + 16 inclinados)")
    print(f"Muestras por conductor: {COIL_COLLISION_CHECK_SAMPLES}")
    print(f"Seccion exterior de carcasa: {2*COIL_CONDUCTOR_RADIUS*1000:.0f} mm radial x {2*COIL_TRANSVERSE_HALF_WIDTH*1000:.0f} mm transversal")
    print(f"Par con menor holgura conservadora: {worst_pair[0]} / {worst_pair[1]}")
    print(f"Holgura superficial conservadora mínima: {worst_gap*1000:.2f} mm")
    print(f"Peor holgura de poste respecto a otra bobina: {worst_support[0]} / {worst_support[1]} = {worst_support_gap*1000:.2f} mm")
    print(f"Holgura mínima exigida: {COIL_MIN_SURFACE_CLEARANCE*1000:.2f} mm")
    print("Holguras de bastidores: " + ", ".join(f"{k}={v*1000:.1f} mm" for k, v in frame_gaps.items()))
    if clashes or support_clashes:
        print("[ERROR] Posibles interferencias detectadas:")
        for a, b, gap in (clashes + support_clashes)[:30]:
            print(f"  - {a} / {b}: holgura conservadora {gap*1000:.2f} mm")
        raise RuntimeError("La verificación previa detectó bobinas o postes demasiado próximos; no se construyó el reactor.")
    print("[OK] No se detectan solapamientos conductor-conductor, poste-conductor ni bastidor-paquete para las secciones elipticas.")
    print("[NOTA] La comprobacion cubre paquetes, postes H y holguras de bastidores; no es una validacion magnetica ni estructural.")
    print("=============================================================\n")
    return worst_gap, worst_pair


def _hatch_gap_phi():
    # La hélice H01 avanza 4 theta por phi; se corta solo donde atraviesa la ventana angular.
    return MAINTENANCE_HATCH_HALF_THETA / POLOIDAL_TURNS_PER_TOROIDAL_TURN + math.radians(1.0)


def _make_helical_points_phi(phi0, phi1, count, phase, minor_radius):
    return [torus_point(phase + POLOIDAL_TURNS_PER_TOROIDAL_TURN * (phi0 + (phi1-phi0)*i/count),
                        phi0 + (phi1-phi0)*i/count, minor_radius) for i in range(count + 1)]


def _make_x_points_theta(theta0, theta1, count, phi_center, direction, is_equatorial, minor_radius):
    pts = []
    for i in range(count + 1):
        theta = theta0 + (theta1-theta0)*i/count
        phi_offset = direction * COIL_X_TILT * (math.sin(theta) if is_equatorial else math.cos(theta))
        pts.append(torus_point(theta, phi_center + phi_offset, minor_radius))
    return pts


def _set_open_mode_hidden(obj):
    obj.hide_set(True)
    obj.hide_render = True


def build_actual_coil_conductors(cols, mats):
    """Crea 4 helicoidales y 16 X; dos segmentos de cruce tienen puente desmontable para servicio."""
    helical_samples = 2048
    for track in range(N_HELICAL_COILS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_HELICAL_COILS
        if track == MAINTENANCE_TRACK:
            gap_half = _hatch_gap_phi()
            gap_start = MAINTENANCE_HATCH_PHI - gap_half
            gap_end = MAINTENANCE_HATCH_PHI + gap_half
            main_count = max(256, int(helical_samples * (2.0*math.pi - 2.0*gap_half) / (2.0*math.pi)))
            bridge_count = max(32, int(helical_samples * (2.0*gap_half) / (2.0*math.pi)))
            points_main = _make_helical_points_phi(gap_end, gap_start + 2.0*math.pi,
                                                   main_count, phase, COIL_CENTERLINE_MINOR)
            points_bridge = _make_helical_points_phi(gap_start, gap_end, bridge_count,
                                                     phase, COIL_CENTERLINE_MINOR)
            obj = create_curve(f"Bobina_Helicoidal_H{track+1:02d}_Tramo_Mantenimiento", points_main,
                               cols["coils"], mats["coil_helix_conductor"], COIL_CORE_RADIAL_HALF_WIDTH,
                               cyclic=False, resolution=1, route_kind="coil", assembly_id=f"H{track+1:02d}")
            create_superconducting_coil_assembly(
                f"Bobina_Helicoidal_H{track+1:02d}_Paquete", points_main, cols["coils"],
                mats["coil_helix_conductor"], cols, mats, f"H{track+1:02d}", cyclic=False)
            bridge_curve = create_curve(f"Puente_Desmontable_Bobina_H{track+1:02d}_Escotilla",
                                        points_bridge, cols["coils"], mats["coil_helix_conductor"],
                                        COIL_CORE_RADIAL_HALF_WIDTH, cyclic=False, resolution=1,
                                        route_kind="coil", assembly_id=f"H{track+1:02d}")
            bridge_parts = create_superconducting_coil_assembly(
                f"Puente_Desmontable_Paquete_Bobina_Helicoidal_H{track+1:02d}_Escotilla",
                points_bridge, cols["coils"], mats["coil_helix_conductor"], cols, mats,
                f"H{track+1:02d}", cyclic=False)
            for bridge in [bridge_curve] + bridge_parts:
                bridge["Pieza_desmontable_para_mantenimiento"] = True
                if MAINTENANCE_HATCH_OPEN:
                    _set_open_mode_hidden(bridge)
            if obj:
                obj["Familia"] = "Helicoidal"
                obj["Vueltas_poloidales_por_toroidal"] = POLOIDAL_TURNS_PER_TOROIDAL_TURN
                obj["Radio_centro_m"] = COIL_CENTERLINE_MINOR
                obj["Radio_conductor_m"] = COIL_CONDUCTOR_RADIUS
                obj["Capa_radial"] = "Capa próxima a la vasija; segmento puente desmontable sobre escotilla T01_C02"
        else:
            points = _make_helical_points_phi(0.0, 2.0*math.pi, helical_samples, phase, COIL_CENTERLINE_MINOR)
            obj = create_curve(f"Bobina_Helicoidal_H{track+1:02d}", points, cols["coils"],
                               mats["coil_helix_conductor"], COIL_CORE_RADIAL_HALF_WIDTH,
                               cyclic=True, resolution=1, route_kind="coil", assembly_id=f"H{track+1:02d}")
            create_superconducting_coil_assembly(
                f"Bobina_Helicoidal_H{track+1:02d}_Paquete", points, cols["coils"],
                mats["coil_helix_conductor"], cols, mats, f"H{track+1:02d}", cyclic=True)
            if obj:
                obj["Familia"] = "Helicoidal"
                obj["Vueltas_poloidales_por_toroidal"] = POLOIDAL_TURNS_PER_TOROIDAL_TURN
                obj["Radio_centro_m"] = COIL_CENTERLINE_MINOR
                obj["Radio_conductor_m"] = COIL_CONDUCTOR_RADIUS
                obj["Capa_radial"] = "Capa próxima a la vasija; separada de la familia X"

    # Las 16 X se mantienen en sus dos carriles habituales; X02B lleva un puente desmontable.
    x_samples = 2048
    for pos in range(8):
        phi_center = pos * (2.0 * math.pi / 8.0)
        is_equatorial = (pos % 2 == 0)
        for coil_index, direction in (("A", 1.0), ("B", -1.0)):
            r_actual = COIL_X_CENTERLINE_MINOR + (
                COIL_X_SEPARATION_R if coil_index == "A" else -COIL_X_SEPARATION_R
            )
            if pos == 1 and coil_index == "B":
                gap_center = MAINTENANCE_HATCH_THETA
                gap_half = MAINTENANCE_HATCH_HALF_THETA + math.radians(4.0)
                gap_start = gap_center - gap_half
                gap_end = gap_center + gap_half
                main_count = max(256, int(x_samples * (2.0*math.pi - 2.0*gap_half) / (2.0*math.pi)))
                bridge_count = max(32, int(x_samples * (2.0*gap_half) / (2.0*math.pi)))
                points_main = _make_x_points_theta(gap_end, gap_start + 2.0*math.pi,
                                                   main_count, phi_center, direction,
                                                   is_equatorial, r_actual)
                points_bridge = _make_x_points_theta(gap_start, gap_end, bridge_count,
                                                     phi_center, direction, is_equatorial, r_actual)
                mat = mats["coil_x_pol_b"]
                obj = create_curve("Bobina_X_X02_B_Tramo_Mantenimiento", points_main, cols["coils"],
                                   mat, COIL_CORE_RADIAL_HALF_WIDTH, cyclic=False, resolution=1, route_kind="coil",
                                   assembly_id="X02_B")
                create_superconducting_coil_assembly(
                    "Bobina_X_X02_B_Paquete", points_main, cols["coils"], mat, cols, mats,
                    "X02_B", cyclic=False)
                bridge_curve = create_curve("Puente_Desmontable_Bobina_X02_B_Escotilla", points_bridge,
                                            cols["coils"], mat, COIL_CORE_RADIAL_HALF_WIDTH,
                                            cyclic=False, resolution=1, route_kind="coil", assembly_id="X02_B")
                bridge_parts = create_superconducting_coil_assembly(
                    "Puente_Desmontable_Paquete_Bobina_X02_B_Escotilla", points_bridge,
                    cols["coils"], mat, cols, mats, "X02_B", cyclic=False)
                for bridge in [bridge_curve] + bridge_parts:
                    bridge["Pieza_desmontable_para_mantenimiento"] = True
                    if MAINTENANCE_HATCH_OPEN:
                        _set_open_mode_hidden(bridge)
                obj["Familia"] = "Superior_Inferior"
                obj["Par_X"] = pos + 1
                obj["Rama"] = coil_index
                obj["Radio_centro_m"] = r_actual
                obj["Inclinacion_maxima_grados"] = math.degrees(COIL_X_TILT)
                obj["Capa_radial"] = "Exterior; tramo de servicio con puente desmontable"
                continue

            points = []
            for j in range(x_samples):
                theta = 2.0 * math.pi * j / x_samples
                if is_equatorial:
                    phi_offset = direction * COIL_X_TILT * math.sin(theta)
                else:
                    phi_offset = direction * COIL_X_TILT * math.cos(theta)
                phi = phi_center + phi_offset
                points.append(torus_point(theta, phi, r_actual))
            if is_equatorial:
                mat = mats["coil_x_eq_a"] if coil_index == "A" else mats["coil_x_eq_b"]
                family = "Ecuatorial"
            else:
                mat = mats["coil_x_pol_a"] if coil_index == "A" else mats["coil_x_pol_b"]
                family = "Superior_Inferior"
            obj = create_curve(f"Bobina_X_{pos+1:02d}_{coil_index}", points, cols["coils"], mat,
                               COIL_CORE_RADIAL_HALF_WIDTH, cyclic=True, resolution=1, route_kind="coil",
                               assembly_id=f"X{pos+1:02d}_{coil_index}")
            create_superconducting_coil_assembly(
                f"Bobina_X_{pos+1:02d}_{coil_index}_Paquete", points, cols["coils"], mat,
                cols, mats, f"X{pos+1:02d}_{coil_index}", cyclic=True)
            if obj:
                obj["Familia"] = family
                obj["Par_X"] = pos + 1
                obj["Rama"] = coil_index
                obj["Radio_centro_m"] = r_actual
                obj["Radio_conductor_m"] = COIL_CONDUCTOR_RADIUS
                obj["Inclinacion_maxima_grados"] = math.degrees(COIL_X_TILT)
                obj["Capa_radial"] = "Exterior; separada de la familia helicoidal"


def x_support_phase_degrees(position_index, coil_index):
    """Fases de apoyo elegidas para que los postes eviten cruzar las helicoidales."""
    if position_index % 2 == 0:  # familia ecuatorial: phi_offset ~ sin(theta)
        return 31.5 if coil_index == "A" else 45.0
    # familia superior/inferior: phi_offset ~ cos(theta)
    return 54.25 if coil_index == "A" else 35.75


def build_coil_supports(cols, mats):
    """Bastidores toroidales y puntales locales para las tres capas superconductoras.

    H queda arriostrada entre la vasija y el aro intermedio. X_B se apoya desde
    el lado interior de ese aro, y X_A desde el lado interior del aro exterior;
    ya no hay postes radiales largos que atraviesen bobinas de otras capas.
    """
    inner_frame_inner = COIL_INNER_FRAME_MINOR - COIL_INNER_FRAME_HALF_THICKNESS
    inner_frame_outer = COIL_INNER_FRAME_MINOR + COIL_INNER_FRAME_HALF_THICKNESS
    outer_frame_inner = COIL_OUTER_FRAME_MINOR - COIL_OUTER_FRAME_HALF_THICKNESS

    inner_frame = create_toroidal_support_shell(
        "Bastidor_Toroidal_Intermedio_Bobinas_H_XB", COIL_INNER_FRAME_MINOR,
        COIL_INNER_FRAME_HALF_THICKNESS, cols["coil_mounts"], mats["support"], "INTERMEDIO",
    )
    outer_frame = create_toroidal_support_shell(
        "Bastidor_Toroidal_Exterior_Bobinas_XA", COIL_OUTER_FRAME_MINOR,
        COIL_OUTER_FRAME_HALF_THICKNESS, cols["coil_mounts"], mats["support"], "EXTERIOR",
    )

    # H1-H4: zapatas/postes cortos desde la vasija a la cara interior de cada hélice.
    for track in range(N_HELICAL_COILS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_HELICAL_COILS
        for seg in range(CASSETTES_PER_TRACK):
            phi = 2.0 * math.pi * (seg + 0.5) / CASSETTES_PER_TRACK
            theta = wrap_angle(phase + POLOIDAL_TURNS_PER_TOROIDAL_TURN * phi)
            if MAINTENANCE_HATCH_OPEN and track == MAINTENANCE_TRACK and seg == MAINTENANCE_SEGMENT:
                continue
            create_coil_support_post(f"H{track+1:02d}_Silla_{seg+1:02d}", theta, phi,
                                     COIL_CENTERLINE_MINOR, COIL_CONDUCTOR_RADIUS,
                                     COIL_HELICAL_SADDLE_SIZE, mats["coil_helix"], cols, mats)

        # Cuatro puntales exteriores H -> cara interior del aro intermedio.
        for k in range(4):
            phi = k * math.pi / 2.0 + math.pi / 8.0
            theta = wrap_angle(phase + POLOIDAL_TURNS_PER_TOROIDAL_TURN * phi)
            if MAINTENANCE_HATCH_OPEN and track == MAINTENANCE_TRACK and angular_distance(phi, MAINTENANCE_HATCH_PHI) < _hatch_gap_phi():
                continue
            create_coil_frame_link(
                f"Puntal_H{track+1:02d}_A_Bastidor_Intermedio_{k+1}", theta, phi,
                COIL_CENTERLINE_MINOR + COIL_CONDUCTOR_RADIUS,
                inner_frame_inner, "INTERMEDIO", f"H{track+1:02d}", cols, mats,
            )

    # X_B: el puntal ocupa exclusivamente el hueco entre H y la capa X_B.
    # X_A: el puntal ocupa exclusivamente el hueco entre la capa X_A y el aro exterior.
    for pos in range(8):
        phi_center = pos * 2.0 * math.pi / 8.0
        is_equatorial = (pos % 2 == 0)
        for coil_index, direction in (("A", 1.0), ("B", -1.0)):
            r_actual = COIL_X_CENTERLINE_MINOR + (
                COIL_X_SEPARATION_R if coil_index == "A" else -COIL_X_SEPARATION_R
            )
            base_theta = math.radians(x_support_phase_degrees(pos, coil_index))
            for sample in range(4):
                theta = base_theta + sample * math.pi / 2.0
                if is_equatorial:
                    phi = phi_center + direction * COIL_X_TILT * math.sin(theta)
                else:
                    phi = phi_center + direction * COIL_X_TILT * math.cos(theta)
                if coil_index == "B":
                    start_minor = inner_frame_outer
                    end_minor = r_actual - COIL_CONDUCTOR_RADIUS
                    frame_id = "INTERMEDIO"
                else:
                    start_minor = r_actual + COIL_CONDUCTOR_RADIUS
                    end_minor = outer_frame_inner
                    frame_id = "EXTERIOR"
                create_coil_frame_link(
                    f"Puntal_X{pos+1:02d}_{coil_index}_A_Bastidor_{sample+1}",
                    wrap_angle(theta), wrap_angle(phi), start_minor, end_minor,
                    frame_id, f"X{pos+1:02d}_{coil_index}", cols, mats,
                )

    create_surface_text("Etiqueta_Soportes_Helicoidales", "BOBINAS HELICOIDALES H1-H4",
                        math.radians(65), math.radians(18), COIL_LABEL_MINOR,
                        cols["labels"], mats["coil_helix"], size=0.018)
    create_surface_text("Etiqueta_Soportes_X", "BOBINAS INCLINADAS X1-X16",
                        math.radians(65), math.radians(202), COIL_LABEL_MINOR,
                        cols["labels"], mats["coil_x"], size=0.018)

# =============================================================================
# 7. PLASMA DE REFERENCIA
# =============================================================================
def create_plasma_reference(collection, mat):
    nphi, ntheta = 384, 128
    verts, faces = [], []
    for i in range(nphi):
        phi = 2.0 * math.pi * i / nphi
        for j in range(ntheta):
            theta = 2.0 * math.pi * j / ntheta
            verts.append(plasma_point(theta, phi))
    for i in range(nphi):
        ni = (i + 1) % nphi
        for j in range(ntheta):
            nj = (j + 1) % ntheta
            a, b = i * ntheta + j, ni * ntheta + j
            c, d = ni * ntheta + nj, i * ntheta + nj
            faces.append((a, b, c, d))
    obj = link_mesh_object("Plasma_Referencia_Helicoidal", verts, faces, collection, mat)
    obj["Nota"] = "Superficie geométrica de referencia; no es equilibrio MHD."
    return obj

# =============================================================================
# 8. CASSETTES, RESPALDOS Y BLANCOS SEGMENTADOS
# =============================================================================
def cassette_phi_bounds(segment):
    center = (segment + 0.5) * TOROIDAL_PITCH
    half = 0.5 * CASSETTE_SPAN
    return center, center - half, center + half


def tile_patch(name, phi0, phi1, th0, th1, phase, collection, mat):
    pm = (phi1 - phi0) * TILE_GAP_FRACTION / 2.0
    tm = (th1 - th0) * TILE_GAP_FRACTION / 2.0
    return create_solid_patch(
        name, phi0 + pm, phi1 - pm, th0 + tm, th1 - tm,
        phase, TILE_INNER_OFFSET, TILE_OUTER_OFFSET,
        collection, mat, n_phi=4, n_theta=3,
    )


def build_cassette(track, segment, phase, cols, mats):
    center, phi_min, phi_max = cassette_phi_bounds(segment)
    base = f"T{track+1:02d}_C{segment+1:02d}"
    # Dos respaldos separados: la ranura central queda realmente abierta.
    wings = (
        (-TARGET_HALF_ANGLE - BACKING_MARGIN_ANGLE, -PUMPING_GAP_HALF, "A"),
        (PUMPING_GAP_HALF, TARGET_HALF_ANGLE + BACKING_MARGIN_ANGLE, "B"),
    )
    for lo, hi, wing in wings:
        create_solid_patch(
            f"{base}_Respaldo_{wing}",
            phi_min - BACKING_PHI_MARGIN, phi_max + BACKING_PHI_MARGIN,
            lo, hi, phase, BACKING_INNER_OFFSET, BACKING_OUTER_OFFSET,
            cols["backing"], mats["backing"], n_phi=10, n_theta=8,
        )
    # 4 columnas toroidales x 2 filas poloidales x 2 alas = 16 teselas/cassette.
    phi_edges = [phi_min + CASSETTE_SPAN * k / TILES_TOROIDAL for k in range(TILES_TOROIDAL + 1)]
    for wing_index, (lo, hi, wing) in enumerate((
        (-TARGET_HALF_ANGLE, -PUMPING_GAP_HALF, "A"),
        (PUMPING_GAP_HALF, TARGET_HALF_ANGLE, "B"),
    )):
        theta_edges = [lo + (hi - lo) * k / TILE_ROWS_PER_WING for k in range(TILE_ROWS_PER_WING + 1)]
        for row_i in range(TILE_ROWS_PER_WING):
            for col_i in range(TILES_TOROIDAL):
                m = mats["tungsten"][(track + segment + wing_index + row_i + col_i) % len(mats["tungsten"])]
                tile_patch(
                    f"{base}_Blanco_{wing}_F{row_i+1}_C{col_i+1}",
                    phi_edges[col_i], phi_edges[col_i + 1],
                    theta_edges[row_i], theta_edges[row_i + 1], phase,
                    cols["targets"], m,
                )
    # Baffles perimetrales elevados; se solapan 0.25 grados con las teselas para
    # no dejar una línea directa de visión detrás del borde.
    for edge, sign in ((-TARGET_HALF_ANGLE, -1), (TARGET_HALF_ANGLE, 1)):
        t0 = edge - BACKING_MARGIN_ANGLE if sign < 0 else edge - BAFFLE_EDGE_OVERLAP
        t1 = edge + BAFFLE_EDGE_OVERLAP if sign < 0 else edge + BACKING_MARGIN_ANGLE
        create_solid_patch(
            f"{base}_Baffle_Perimetral_{'Inf' if sign < 0 else 'Sup'}",
            phi_min, phi_max, t0, t1, phase,
            BAFFLE_INNER_OFFSET, BAFFLE_OUTER_OFFSET,
            cols["baffles"], mats["baffle"], n_phi=8, n_theta=3,
        )
    # Cierres toroidales breves de las alas. No invaden el canal de servicios posterior.
    for edge_i, (p0, p1) in enumerate(((phi_min - END_FLANGE_PHI, phi_min), (phi_max, phi_max + END_FLANGE_PHI))):
        for wing, lo, hi in (
            ("A", -TARGET_HALF_ANGLE - BACKING_MARGIN_ANGLE, -PUMPING_GAP_HALF),
            ("B", PUMPING_GAP_HALF, TARGET_HALF_ANGLE + BACKING_MARGIN_ANGLE),
        ):
            create_solid_patch(
                f"{base}_Cierre_{edge_i+1}_{wing}", p0, p1, lo, hi, phase,
                BAFFLE_INNER_OFFSET, BAFFLE_OUTER_OFFSET,
                cols["baffles"], mats["baffle"], n_phi=2, n_theta=5,
            )
    return phi_min, phi_max

# =============================================================================
# 9. PLENUM HELICOIDAL: SUELO CON ABERTURAS Y PAREDES LATERALES
# =============================================================================
def build_pumping_plenum(track, phase, cols, mats):
    """Suelo y paredes del plenum; T01_C02 incorpora insertos desmontables de servicio."""
    maintenance_for_this_track = (track == MAINTENANCE_TRACK)
    maintenance_hole_floor = []
    if maintenance_for_this_track:
        maintenance_hole_floor.append((MAINTENANCE_HATCH_PHI, 0.0,
                                       MAINTENANCE_HATCH_HALF_PHI + 0.007,
                                       PUMPING_GAP_HALF + math.radians(0.6)))
    create_cyclic_solid_ribbon(
        f"T{track+1:02d}_Suelo_Pleno_con_Boca",
        phase, -PUMPING_GAP_HALF, PUMPING_GAP_HALF,
        PLENUM_FLOOR_INNER_OFFSET, PLENUM_FLOOR_OUTER_OFFSET,
        cols["baffles"], mats["plenum"], n_phi=PLENUM_N_PHI, n_theta=PLENUM_N_THETA,
        hole=(pump_port_phi(track, phase), 0.0, PLENUM_HOLE_HALF_PHI, PLENUM_HOLE_HALF_THETA),
        holes=maintenance_hole_floor,
    )
    for side, center in (("A", -PUMPING_GAP_HALF), ("B", PUMPING_GAP_HALF)):
        half = SLOT_WALL_THICKNESS_ANGLE / 2.0
        wall_holes = []
        if maintenance_for_this_track:
            wall_holes.append((MAINTENANCE_HATCH_PHI, center,
                               MAINTENANCE_HATCH_HALF_PHI + 0.007,
                               math.radians(0.60)))
        create_cyclic_solid_ribbon(
            f"T{track+1:02d}_Pared_Lateral_Pleno_{side}",
            phase, center - half, center + half,
            SLOT_WALL_INNER_OFFSET, SLOT_WALL_OUTER_OFFSET,
            cols["baffles"], mats["baffle"], n_phi=PLENUM_N_PHI, n_theta=2,
            holes=wall_holes,
        )
    if maintenance_for_this_track:
        # Insertos que cierran el tramo de plenum mientras el reactor está en configuración normal.
        floor_insert = create_solid_patch(
            "Inserto_Desmontable_Suelo_Pleno_T01_C02",
            MAINTENANCE_HATCH_PHI - MAINTENANCE_HATCH_HALF_PHI - 0.007,
            MAINTENANCE_HATCH_PHI + MAINTENANCE_HATCH_HALF_PHI + 0.007,
            -PUMPING_GAP_HALF - math.radians(0.6), PUMPING_GAP_HALF + math.radians(0.6),
            phase, PLENUM_FLOOR_INNER_OFFSET, PLENUM_FLOOR_OUTER_OFFSET,
            cols["maintenance"], mats["maintenance"], n_phi=14, n_theta=4,
        )
        floor_insert["Desmontable"] = True
        if MAINTENANCE_HATCH_OPEN:
            floor_insert.hide_set(True)
            floor_insert.hide_render = True
        for side, center in (("A", -PUMPING_GAP_HALF), ("B", PUMPING_GAP_HALF)):
            insert = create_solid_patch(
                f"Inserto_Desmontable_Pared_Pleno_{side}_T01_C02",
                MAINTENANCE_HATCH_PHI - MAINTENANCE_HATCH_HALF_PHI - 0.007,
                MAINTENANCE_HATCH_PHI + MAINTENANCE_HATCH_HALF_PHI + 0.007,
                center - half, center + half, phase,
                SLOT_WALL_INNER_OFFSET, SLOT_WALL_OUTER_OFFSET,
                cols["maintenance"], mats["maintenance"], n_phi=14, n_theta=2,
            )
            insert["Desmontable"] = True
            if MAINTENANCE_HATCH_OPEN:
                insert.hide_set(True)
                insert.hide_render = True

# =============================================================================
# 10. REFRIGERACIÓN DE LOS CASSETTES Y COLECTORES
# =============================================================================
def serpentine_points(phase, phi_min, phi_max, theta_lo, theta_hi):
    """Serpentin de cuatro pasadas; ambos extremos quedan en el extremo phi_min."""
    fracs = (0.15, 0.38, 0.62, 0.85)
    rows = [theta_lo + (theta_hi - theta_lo) * f for f in fracs]
    # Dejar 4 mm geométricos de margen respecto a los extremos del cassette.
    p0 = phi_min + 0.0040
    p1 = phi_max - 0.0040
    pts = []
    for row_i, thoff in enumerate(rows):
        indices = range(CHANNEL_SAMPLES_PER_PASS) if row_i % 2 == 0 else range(CHANNEL_SAMPLES_PER_PASS - 1, -1, -1)
        pass_pts = []
        for idx in indices:
            f = idx / (CHANNEL_SAMPLES_PER_PASS - 1)
            phi = p0 + (p1 - p0) * f
            pass_pts.append(make_track_point(phi, thoff, phase, CHANNEL_CENTER_OFFSET))
        if pts:
            # Conector curvo en el extremo de la pasada; evita saltos entre filas.
            prev = pts[-1]
            nxt = pass_pts[0]
            phi_turn = p1 if row_i % 2 == 1 else p0
            last_off = rows[row_i - 1]
            for k in range(1, 25):
                f = k / 25.0
                th = last_off + (thoff - last_off) * f
                pts.append(make_track_point(phi_turn, th, phase, CHANNEL_CENTER_OFFSET))
        pts.extend(pass_pts)
    return pts, rows, p0, p1


def create_cooling_branch(name, phase, phi_start, theta_source, theta_header, stage_delta,
                          target_minor_radius, collection, mat):
    """Conecta el serpentín con un carril radial exclusivo de su colector.

    El ascenso inicial parte del extremo real del serpentín; después la ruta usa
    un radio menor absoluto fijo, por lo que la ondulación del plasma no acerca
    los ramales a los colectores vecinos.
    """
    theta0 = ribbon_theta(phi_start, theta_source, phase)
    start_minor = plasma_minor_radius(theta0, phi_start) + CHANNEL_CENTER_OFFSET
    points = [torus_point(theta0, phi_start, start_minor)]
    # Elevación radial desde el canal al carril asignado; densificación para chequeos.
    for k in range(1, 25):
        f = k / 24.0
        minor = start_minor + (target_minor_radius - start_minor) * f
        points.append(torus_point(theta0, phi_start, minor))

    phi_stage = phi_start + stage_delta
    # Pequeño desplazamiento toroidal con theta relativo constante.
    for k in range(1, 101):
        f = k / 100.0
        phi = phi_start + (phi_stage - phi_start) * f
        theta = ribbon_theta(phi, theta_source, phase)
        points.append(torus_point(theta, phi, target_minor_radius))

    # Transición poloidal en el radio propio de este ramal.
    for k in range(1, 201):
        f = k / 200.0
        theta_off = theta_source + (theta_header - theta_source) * f
        theta = ribbon_theta(phi_stage, theta_off, phase)
        points.append(torus_point(theta, phi_stage, target_minor_radius))

    return create_curve(name, points, collection, mat, BRANCH_RADIUS, cyclic=False, resolution=1,
                        route_kind="coolant")


def build_cassette_cooling(track, segment, phase, phi_min, phi_max, cols, mats):
    """Un serpentín continuo por ala, integrado en el reverso del respaldo.

    Las conexiones se enrutan por cuatro carriles separados en el hueco toroidal
    entre cassettes, evitando duplicar tubos sobre el mismo serpentín.
    """
    base = f"T{track+1:02d}_C{segment+1:02d}"
    wings = (
        {
            "name": "A", "lo": -TARGET_HALF_ANGLE, "hi": -PUMPING_GAP_HALF,
            "supply_off": -HEADER_THETA_FAR, "supply_lane": -0.0150,
            "return_off": -HEADER_THETA_NEAR, "return_lane": -0.0150,
        },
        {
            "name": "B", "lo": PUMPING_GAP_HALF, "hi": TARGET_HALF_ANGLE,
            "supply_off": HEADER_THETA_NEAR, "supply_lane": 0.0150,
            "return_off": HEADER_THETA_FAR, "return_lane": 0.0150,
        },
    )
    for wing in wings:
        pts, rows, p0, p1 = serpentine_points(
            phase, phi_min, phi_max, wing["lo"], wing["hi"]
        )
        # Un solo canal físico por ala: entrada y salida son sus dos extremos.
        create_curve(
            f"{base}_Serpentin_{wing['name']}", pts, cols["cooling"], mats["channel"],
            CHANNEL_RADIUS, cyclic=False, resolution=1, route_kind="coolant",
        )
        create_cooling_branch(
            f"{base}_Ramal_{wing['name']}_Ida", phase, p0, rows[0],
            wing["supply_off"], wing["supply_lane"],
            coolant_header_radius(f"{wing['name']}_Ida"), cols["cooling"], mats["supply"],
        )
        create_cooling_branch(
            f"{base}_Ramal_{wing['name']}_Retorno", phase, p0, rows[-1],
            wing["return_off"], wing["return_lane"],
            coolant_header_radius(f"{wing['name']}_Retorno"), cols["cooling"], mats["return"],
        )


def build_helical_headers(track, phase, cols, mats):
    # Cada ala tiene su colector de ida y retorno. Los cuatro trazados por banda están
    # separados angularmente para dejar libre el plenum central y los apoyos.
    definitions = (
        ("A_Ida", -HEADER_THETA_FAR, mats["supply"]),
        ("A_Retorno", -HEADER_THETA_NEAR, mats["return"]),
        ("B_Ida", HEADER_THETA_NEAR, mats["supply"]),
        ("B_Retorno", HEADER_THETA_FAR, mats["return"]),
    )
    for label, theta_off, mat in definitions:
        minor_radius = COOLANT_COLLECTOR_RADII[label]
        # 8192 nodos reducen la separación entre vértices de la polilínea tubular.
        pts = [torus_point(ribbon_theta(2.0 * math.pi * i / 8192, theta_off, phase),
                           2.0 * math.pi * i / 8192, minor_radius)
               for i in range(8192)]
        create_curve(
            f"T{track+1:02d}_Colector_{label}", pts, cols["headers"], mat,
            HEADER_RADIUS, cyclic=True, resolution=1, route_kind="coolant",
        )

# =============================================================================
# 11. TUBOS DE BOMBEO, PASOS Y SOPORTES
# =============================================================================
def build_pump_duct(track, phase, cols, mats):
    """Conducto de vacío guiado por encima de los ramales de agua.

    Sale por el centro de la ranura, sube radialmente al carril de vacío y solo
    entonces se desplaza 18° hacia el pasamuros. Así no cruza los serpentines.
    """
    phi = pump_port_phi(track, phase)
    points, radii = [], []
    theta0 = ribbon_theta(phi, 0.0, phase)
    start_minor = plasma_minor_radius(theta0, phi) + PLENUM_FLOOR_OUTER_OFFSET
    taper_ratio = PUMP_DUCT_MOUTH_RADIUS / PUMP_DUCT_RADIUS
    # Boca ensanchada en el suelo del plenum que se reduce al subir al carril.
    for k in range(81):
        f = k / 80.0
        minor = start_minor + (PUMP_SERVICE_MINOR_RADIUS - start_minor) * f
        points.append(torus_point(theta0, phi, minor))
        radii.append(taper_ratio - (taper_ratio - 1.0) * f)
    # Desplazamiento poloidal por el carril exterior, fuera de los circuitos de agua.
    for k in range(1, 201):
        f = k / 200.0
        theta_off = PUMP_PORT_THETA_OFFSET * f
        theta = ribbon_theta(phi, theta_off, phase)
        points.append(torus_point(theta, phi, PUMP_SERVICE_MINOR_RADIUS))
        radii.append(1.0)
    # Paso por la pared y extensión exterior hacia la bomba.
    theta_port = ribbon_theta(phi, PUMP_PORT_THETA_OFFSET, phase)
    for k in range(1, 101):
        f = k / 100.0
        minor = PUMP_SERVICE_MINOR_RADIUS + (PUMP_STUB_MINOR - PUMP_SERVICE_MINOR_RADIUS) * f
        points.append(torus_point(theta_port, phi, minor))
        radii.append(1.0)
    create_curve(
        f"T{track+1:02d}_Conducto_Extraccion_Pleno", points, cols["feedthroughs"],
        mats["pump"], PUMP_DUCT_RADIUS, cyclic=False, point_radii=radii, resolution=1, route_kind="vacuum",
    )


def create_supports(track, segment, phase, phi_center, cols, mats):
    base = f"T{track+1:02d}_C{segment+1:02d}"
    for idx, theta_off in enumerate((-SUPPORT_THETA_OFFSET, SUPPORT_THETA_OFFSET), 1):
        phi = phi_center
        theta = ribbon_theta(phi, theta_off, phase)
        p0 = surface_point(theta, phi, BACKING_OUTER_OFFSET - 0.0005)
        p1 = torus_point(theta, phi, SUPPORT_WALL_CENTER_MINOR)
        create_curve(
            f"{base}_Puntal_{idx}", [p0, p1], cols["supports"], mats["support"],
            SUPPORT_RADIUS, cyclic=False, resolution=1,
        )
        create_uv_sphere(
            f"{base}_Anclaje_{idx}", p1, SUPPORT_ANCHOR_RADIUS,
            cols["supports"], mats["support"], segments=12, rings=8,
        )
        # Zapata sobre el reverso del respaldo, en la zona del hombro.
        pad_half_phi = 0.014
        pad_half_th = math.radians(1.3)
        create_solid_patch(
            f"{base}_Zapata_{idx}", phi - pad_half_phi, phi + pad_half_phi,
            theta_off - pad_half_th, theta_off + pad_half_th, phase,
            BACKING_OUTER_OFFSET - 0.003, BACKING_OUTER_OFFSET + 0.002,
            cols["supports"], mats["support"], n_phi=4, n_theta=3,
        )

# =============================================================================
# 12. COMPROBACIÓN INTEGRAL DE TUBERÍAS CONTRA TUBERÍAS Y BOBINAS
# =============================================================================
def densify_polyline(points, scales, base_radius, max_step=PIPE_CHECK_MAX_STEP):
    """Densifica una polilínea e interpola el radio variable de sus puntos."""
    if len(points) < 2:
        pts = [Vector(p).copy() for p in points]
        radii = [base_radius * s for s in scales]
        return pts, radii
    out = [Vector(points[0]).copy()]
    out_radii = [base_radius * scales[0]]
    for idx, (a, b) in enumerate(zip(points[:-1], points[1:])):
        a, b = Vector(a), Vector(b)
        scale_a, scale_b = scales[idx], scales[idx + 1]
        length = (b - a).length
        n = max(1, int(math.ceil(length / max_step)))
        for k in range(1, n + 1):
            f = k / n
            out.append(a.lerp(b, f))
            scale = scale_a + (scale_b - scale_a) * f
            out_radii.append(base_radius * scale)
    return out, out_radii


def intended_pipe_connection(a, b):
    """Excluye solo pares que deben tocarse para formar un circuito conectado."""
    # Ramal con su serpentín del mismo cassette y de la misma ala.
    if ("_Serpentin_" in a and "_Ramal_" in b) or ("_Serpentin_" in b and "_Ramal_" in a):
        serp, branch = (a, b) if "_Serpentin_" in a else (b, a)
        serp_base, wing = serp.split("_Serpentin_", 1)
        branch_base, role = branch.split("_Ramal_", 1)
        return serp_base == branch_base and wing in role
    # Ramal al colector al que está asignado.
    if ("_Ramal_" in a and "_Colector_" in b) or ("_Ramal_" in b and "_Colector_" in a):
        branch, header = (a, b) if "_Ramal_" in a else (b, a)
        branch_track = branch.split("_C", 1)[0]
        role = branch.split("_Ramal_", 1)[1]
        return header == f"{branch_track}_Colector_{role}"
    # Pasamuros de agua con su colector homónimo.
    if a.startswith("Paso_") or b.startswith("Paso_"):
        port, other = (a, b) if a.startswith("Paso_") else (b, a)
        parts = port[len("Paso_"):].split("_", 1)
        if len(parts) == 2 and "_Colector_" in other:
            return other == f"{parts[0]}_Colector_{parts[1]}"
    # Conducto de extracción y enlace inmediato al módulo de vacío.
    if ("_Conducto_Extraccion_Pleno" in a and "_Enlace_Bomba_Vacio" in b) or (
            "_Conducto_Extraccion_Pleno" in b and "_Enlace_Bomba_Vacio" in a):
        return a[:3] == b[:3]
    # Conexiones intencionadas de la precámara de combustible.
    if {a, b} == {"Linea_Evacuacion_Precamara", "Conexion_Evacuacion_Precamara_A_Colector_Vacio"}:
        return True
    if {a, b} == {"Linea_Llenado_Combustible_Precamara", "Conexion_Valvula_Llenado_A_Botella"}:
        return True
    # El lead HV se conecta intencionadamente a su electrodo correspondiente.
    for pole in ("Anodo", "Catodo"):
        if {a, b} == {f"Electrodo_Helicoidal_{pole}_Ensayo", f"Cable_HV_{pole}_A_Ranura_Independiente"}:
            return True
    # Las cuatro bombas turbomoleculares descargan en el colector común; este
    # colector conecta después con la bomba primaria de respaldo.
    if "Colector_Exhaust_Vacio" in a or "Colector_Exhaust_Vacio" in b:
        other = b if "Colector_Exhaust_Vacio" in a else a
        if ("_Escape_Bomba" in other or other == "Linea_Principal_Vacio"
                or other == "Conexion_Evacuacion_Precamara_A_Colector_Vacio"):
            return True
    return False


def path_bbox(points, radius, padding=PIPE_CHECK_MAX_STEP):
    xs = [p.x for p in points]; ys = [p.y for p in points]; zs = [p.z for p in points]
    pad = radius + padding
    return ((min(xs)-pad, min(ys)-pad, min(zs)-pad),
            (max(xs)+pad, max(ys)+pad, max(zs)+pad))


def boxes_disjoint(a, b):
    lo_a, hi_a = a; lo_b, hi_b = b
    return any(hi_a[i] < lo_b[i] or hi_b[i] < lo_a[i] for i in range(3))


def check_all_pipe_clearances():
    """Comprueba todas las polilíneas de agua/vacío frente a bobinas y otras tuberías.

    Usa las mismas coordenadas y radios pasados a create_curve. La holgura que
    informa resta el radio de ambos tubos y media longitud máxima de segmento de
    muestreo por ruta; si la cota inferior no supera 1 mm, se marca para revisión.
    """
    if not CLEARANCE_PATHS:
        print("[AVISO] No hay rutas registradas para la comprobación integral.")
        return

    routes = []
    for spec in CLEARANCE_PATHS:
        pts, radii = densify_polyline(spec["points"], spec["scales"], spec["radius"])
        max_step = max(((b-a).length for a, b in zip(pts[:-1], pts[1:])), default=0.0)
        max_radius = max(radii) if radii else spec["radius"]
        tree = KDTree(len(pts))
        for idx, p in enumerate(pts):
            tree.insert(p, idx)
        tree.balance()
        routes.append({
            "name": spec["name"], "kind": spec["kind"], "radius": spec["radius"],
            "assembly_id": spec.get("assembly_id"),
            "points": pts, "radii": radii, "tree": tree, "max_step": max_step,
            "bbox": path_bbox(pts, max_radius),
            "variable_radius": (max(radii) - min(radii) > 1e-9) if radii else False,
        })

    coils = [r for r in routes if r["kind"] == "coil"]
    tubes = [r for r in routes if r["kind"] != "coil"]
    failures = []
    best_coil_gap = (float("inf"), "", "")
    best_tube_gap = (float("inf"), "", "")

    def conservative_gap(a, b):
        """Distancia inferior entre superficies con radios interpolados por muestra."""
        short, long = (a, b) if len(a["points"]) <= len(b["points"]) else (b, a)
        min_surface = float("inf")
        def scan(src, dst):
            local_min = float("inf")
            for i, p in enumerate(src["points"]):
                _co, j, dist = dst["tree"].find(p)
                gap_here = dist - src["radii"][i] - dst["radii"][j]
                if gap_here < local_min:
                    local_min = gap_here
            return local_min
        min_surface = scan(short, long)
        # Cuando el radio cambia (boca ensanchada de vacío), revisar también desde el otro lado.
        if a["variable_radius"] or b["variable_radius"]:
            min_surface = min(min_surface, scan(long, short))
        return min_surface - 0.5 * (a["max_step"] + b["max_step"])

    # Cada tubo frente a las 20 bobinas: el margen se calcula de forma conservadora.
    for tube in tubes:
        for coil in coils:
            if tube.get("assembly_id") and tube.get("assembly_id") == coil.get("assembly_id"):
                continue  # canal/boquilla interna de su propia bobina: contacto deliberado
            if boxes_disjoint(tube["bbox"], coil["bbox"]):
                continue
            gap = conservative_gap(tube, coil)
            if gap < best_coil_gap[0]:
                best_coil_gap = (gap, tube["name"], coil["name"])
            if gap < PIPE_CLEARANCE_MIN:
                failures.append((gap, tube["name"], coil["name"], "tubo-bobina"))

    # Tubo contra tubo; se omiten solo uniones explícitas del mismo circuito.
    for i, a in enumerate(tubes):
        for b in tubes[i+1:]:
            if (a["kind"] == "cryo_internal" and b["kind"] == "cryo_internal"
                    and a.get("assembly_id") and a.get("assembly_id") == b.get("assembly_id")):
                continue  # ramificaciones y tomas del circuito interno de una misma bobina
            if intended_pipe_connection(a["name"], b["name"]):
                continue
            if boxes_disjoint(a["bbox"], b["bbox"]):
                continue
            gap = conservative_gap(a, b)
            if gap < best_tube_gap[0]:
                best_tube_gap = (gap, a["name"], b["name"])
            if gap < PIPE_CLEARANCE_MIN:
                failures.append((gap, a["name"], b["name"], "tubo-tubo"))

    print("\n========== COMPROBACIÓN DE TUBERÍAS V14 ==========")
    print(f"Rutas examinadas: {len(tubes)} tubos/ramales (incluye helio interno) y {len(coils)} conductores de bobina")
    if best_coil_gap[1]:
        print(f"Menor holgura conservadora tubo-bobina: {best_coil_gap[0]*1000:.2f} mm | {best_coil_gap[1]} / {best_coil_gap[2]}")
    if best_tube_gap[1]:
        print(f"Menor holgura conservadora tubo-tubo: {best_tube_gap[0]*1000:.2f} mm | {best_tube_gap[1]} / {best_tube_gap[2]}")
    print(f"Holgura mínima exigida: {PIPE_CLEARANCE_MIN*1000:.1f} mm")
    if failures:
        print(f"[ERROR] Se detectaron {len(failures)} posibles interferencias:")
        for gap, a, b, kind in sorted(failures)[:40]:
            print(f"  - {kind}: {a} <> {b}: cota inferior {gap*1000:.2f} mm")
        raise RuntimeError("La comprobación integral detectó rutas con holgura insuficiente; revisar antes de usar el modelo.")
    print("[OK] No se detectan solapamientos de rutas de agua/vacío con bobinas ni entre tubos, con la cota conservadora indicada.")
    print("[NOTA] Es una comprobación geométrica de las rutas de curvas; no comprueba todos los sólidos, soportes, soldaduras o tolerancias de fabricación.")
    print("================================================\n")


# =============================================================================
# 12A. PRIMERA PARED, BLANKET CONCEPTUAL Y BLINDAJE NEUTRÓNICO
# =============================================================================
def panel_intersects_any_port(phi0, phi1, theta_off0, theta_off1, phase, ports):
    """Prueba rápida de intersección angular para no colocar paneles sobre puertos."""
    for port in ports:
        # Incluye también la ventana: los módulos interiores no deben tapar su línea de visión.
        # Muestreo mas denso que el tamaño de los puertos estrechos (combustible/ventana)
        # para reducir el riesgo de que una penetración caiga entre dos muestras.
        for i in range(33):
            phi = phi0 + (phi1-phi0)*i/32.0
            if angular_distance(phi, port["phi"]) > port["half_phi"] + 0.002:
                continue
            for j in range(17):
                thoff = theta_off0 + (theta_off1-theta_off0)*j/16.0
                theta = wrap_angle(ribbon_theta(phi, thoff, phase))
                if angular_distance(theta, port["theta"]) <= port["half_theta"] + 0.010:
                    return True
    return False


def panel_bands_avoiding_collectors():
    """Devuelve bandas poloidales para paneles, dejando corredores para colectores.

    Los corredores son continuos alrededor del toro y están centrados en las cuatro
    pistas de colector que cruzan cada hueco entre bandas de divertor. El ancho angular
    incluye el radio del tubo y un margen de servicio; no se colocan paneles dentro.
    """
    lanes = sorted(COLLECTOR_SERVICE_OFFSETS)
    half = COLLECTOR_SERVICE_HALF_ANGLE
    bands = []
    cursor = GAP_PANEL_THETA_LO
    for center in lanes:
        lo = center - half
        hi = center + half
        if lo > cursor + 1e-8:
            bands.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < GAP_PANEL_THETA_HI - 1e-8:
        bands.append((cursor, GAP_PANEL_THETA_HI))
    return [(lo, hi) for lo, hi in bands if hi - lo > math.radians(0.25)]


def collector_service_corridors_are_clear():
    """Verifica que los corredores de colectores no cortan el área de paneles."""
    bands = panel_bands_avoiding_collectors()
    for lane in COLLECTOR_SERVICE_OFFSETS:
        lane_lo = lane - COLLECTOR_SERVICE_HALF_ANGLE
        lane_hi = lane + COLLECTOR_SERVICE_HALF_ANGLE
        if any(max(lane_lo, lo) < min(lane_hi, hi) - 1e-8 for lo, hi in bands):
            return False
    return bool(bands)


def build_first_wall_blanket_shield(cols, mats, ports):
    """Capas segmentadas en las cuatro bandas libres entre las pistas de divertor.

    Primera pared, blanket y blindaje se segmentan alrededor de corredores de
    colectores. Cada corredor queda libre en toda la vuelta toroidal; los paneles
    se reparten entre los corredores y conservan un tamaño modular razonable.
    El mismo keep-out angular se aplica a cada submódulo para no retirar módulos
    completos cuando solo una pequeña parte coincide con una penetración. Las bandas son la
    superficie disponible entre las huellas de los divertores; no se afirma que
    cubran la zona ocupada por los propios cassettes ni que sean un blanket D-T
    o blindaje neutrónico dimensionado.
    """
    tile_count = blanket_count = shield_count = skipped_fw = skipped_back = 0
    panel_bands = panel_bands_avoiding_collectors()
    if not collector_service_corridors_are_clear():
        raise RuntimeError("Los corredores de colectores se solapan con los paneles.")

    for track in range(N_TRACKS):
        phase = THETA_PHASE + 2.0 * math.pi * track / N_TRACKS
        for segment in range(CASSETTES_PER_TRACK):
            phi_seg0 = segment * TOROIDAL_PITCH + FW_PANEL_PHI_MARGIN
            phi_seg1 = (segment + 1) * TOROIDAL_PITCH - FW_PANEL_PHI_MARGIN
            span_phi = phi_seg1 - phi_seg0

            # Primera pared: una banda modular por cada intervalo libre entre colectores.
            for row_i, (th0, th1) in enumerate(panel_bands):
                for col_i in range(FW_TILES_TOROIDAL):
                    p0 = phi_seg0 + span_phi * col_i / FW_TILES_TOROIDAL
                    p1 = phi_seg0 + span_phi * (col_i + 1) / FW_TILES_TOROIDAL
                    pm = (p1 - p0) * FW_TILE_GAP_FRACTION * 0.5
                    tm = (th1 - th0) * FW_TILE_GAP_FRACTION * 0.5
                    a0, a1 = p0 + pm, p1 - pm
                    b0, b1 = th0 + tm, th1 - tm
                    if panel_intersects_any_port(a0, a1, b0, b1, phase, ports):
                        skipped_fw += 1
                        continue
                    tile = create_solid_patch(
                        f"FW_T{track+1:02d}_S{segment+1:02d}_R{row_i+1}_C{col_i+1}",
                        a0, a1, b0, b1, phase,
                        FIRST_WALL_INNER_OFFSET, FIRST_WALL_OUTER_OFFSET,
                        cols["layers"], mats["first_wall"], n_phi=5, n_theta=3,
                    )
                    tile["Subsystem"] = "Primera pared segmentada"
                    tile["Desmontable"] = True
                    tile["Material_conceptual"] = "Tungsteno/aleación; grado y unión no especificados"
                    tile["Banda_poloidal_offset_grados"] = [math.degrees(th0), math.degrees(th1)]
                    tile["Corredores_colectores_excluidos"] = True
                    tile_count += 1

            # Blanket y blindaje siguen exactamente los mismos corredores despejados.
            for row_i, (th0, th1) in enumerate(panel_bands):
                for col_i in range(BACKING_TILES_TOROIDAL):
                    p0 = phi_seg0 + span_phi * col_i / BACKING_TILES_TOROIDAL
                    p1 = phi_seg0 + span_phi * (col_i + 1) / BACKING_TILES_TOROIDAL
                    pm = (p1 - p0) * BACKING_TILE_GAP_FRACTION * 0.5
                    tm = (th1 - th0) * BACKING_TILE_GAP_FRACTION * 0.5
                    a0, a1 = p0 + pm, p1 - pm
                    b0, b1 = th0 + tm, th1 - tm
                    if panel_intersects_any_port(a0, a1, b0, b1, phase, ports):
                        skipped_back += 1
                        continue
                    blanket = create_solid_patch(
                        f"Blanket_T{track+1:02d}_S{segment+1:02d}_R{row_i+1}_C{col_i+1}",
                        a0, a1, b0, b1, phase,
                        BLANKET_INNER_OFFSET, BLANKET_OUTER_OFFSET,
                        cols["layers"], mats["blanket"], n_phi=8, n_theta=4,
                    )
                    blanket["Subsystem"] = "Blanket modular conceptual D-T"
                    blanket["Desmontable"] = True
                    blanket["Nota"] = "Sin breeder de litio dimensionado ni validacion neutrónica/termica"
                    blanket["Banda_poloidal_offset_grados"] = [math.degrees(th0), math.degrees(th1)]
                    blanket["Corredores_colectores_excluidos"] = True
                    blanket_count += 1

                    shield = create_solid_patch(
                        f"Blindaje_Neutronico_T{track+1:02d}_S{segment+1:02d}_R{row_i+1}_C{col_i+1}",
                        a0, a1, b0, b1, phase,
                        NEUTRON_SHIELD_INNER_OFFSET, NEUTRON_SHIELD_OUTER_OFFSET,
                        cols["layers"], mats["neutron_shield"], n_phi=8, n_theta=4,
                    )
                    shield["Subsystem"] = "Blindaje neutrónico modular conceptual"
                    shield["Desmontable"] = True
                    shield["Espesor_nominal_m"] = NEUTRON_SHIELD_OUTER_OFFSET - NEUTRON_SHIELD_INNER_OFFSET
                    shield["Banda_poloidal_offset_grados"] = [math.degrees(th0), math.degrees(th1)]
                    shield["Corredores_colectores_excluidos"] = True
                    shield_count += 1

            # Dos anclajes por sector solo cuando sus puntos de salida no coinciden con puertos.
            for side, ph in enumerate((phi_seg0 + 0.22 * span_phi, phi_seg0 + 0.78 * span_phi), 1):
                # Anclajes centrados en la banda media, lejos de los carriles de colector.
                thoff = 0.5 * (GAP_PANEL_THETA_LO + GAP_PANEL_THETA_HI)
                th = ribbon_theta(ph, thoff, phase)
                near_port = any(
                    angular_distance(ph, port["phi"]) <= port["half_phi"] + 0.004 and
                    angular_distance(th, port["theta"]) <= port["half_theta"] + math.radians(1.0)
                    for port in ports
                )
                if near_port:
                    continue
                p_start = surface_point(th, ph, NEUTRON_SHIELD_OUTER_OFFSET)
                p_end = torus_point(th, ph, R_MINOR_CHAMBER - 0.0015)
                create_curve(f"Anclaje_Blindaje_T{track+1:02d}_S{segment+1:02d}_{side}",
                             [p_start, p_end], cols["layers"], mats["support"], 0.0012, cyclic=False)

    fw_nominal = N_TRACKS * CASSETTES_PER_TRACK * FW_TILES_TOROIDAL * len(panel_bands)
    back_nominal = N_TRACKS * CASSETTES_PER_TRACK * BACKING_TILES_TOROIDAL * len(panel_bands)
    print(f"[V15] Primera pared: {tile_count}/{fw_nominal} submódulos; omitidos por portales: {skipped_fw}.")
    print(f"[V14] Blanket: {blanket_count}/{back_nominal}; blindaje: {shield_count}/{back_nominal}; "
          f"submódulos omitidos por portales/accesos: {skipped_back}.")
    print("[V14] Corredores de colectores libres en toda la vuelta: " + ", ".join(
        f"{v:.2f}°" for v in map(math.degrees, COLLECTOR_SERVICE_OFFSETS)))
    print("[V14] Cobertura en bandas libres entre los divertores, con corredores de colectores reservados; no es cobertura integral de 360°.")


# =============================================================================
# 12B. GEOMETRÍA AUXILIAR DEL INYECTOR Y ELECTRODOS DE ENSAYO HV
# =============================================================================
def create_oriented_hollow_tube(name, theta, phi, minor_start, minor_end,
                                outer_radius, inner_radius, collection, mat, segments=48):
    """Tubo cilíndrico hueco alineado con la normal local del toro, con ambos extremos anulares."""
    e1, e2, normal = torus_basis(theta, phi)  # normal documentada para mantener el sistema de ejes local
    verts = []
    # Orden: extremo inicio exterior/interior, extremo final exterior/interior.
    for minor, radius in ((minor_start, outer_radius), (minor_start, inner_radius),
                          (minor_end, outer_radius), (minor_end, inner_radius)):
        c = torus_point(theta, phi, minor)
        for i in range(segments):
            a = 2.0*math.pi*i/segments
            verts.append(c + e1*(radius*math.cos(a)) + e2*(radius*math.sin(a)))
    def ix(layer, i): return layer*segments + (i % segments)
    faces=[]
    for i in range(segments):
        j=(i+1)%segments
        faces.append((ix(0,i),ix(0,j),ix(2,j),ix(2,i)))
        faces.append((ix(1,j),ix(1,i),ix(3,i),ix(3,j)))
        faces.append((ix(0,j),ix(0,i),ix(1,i),ix(1,j)))
        faces.append((ix(2,i),ix(2,j),ix(3,j),ix(3,i)))
    obj=link_mesh_object(name,verts,faces,collection,mat,smooth=True)
    obj["Hueco_interior_m"] = inner_radius
    obj["Tipo"] = "Carcasa/tubo conceptual; requiere cálculo de presión y material"
    return obj


def _tube_points_at_minor(theta, phi, r0, r1, count=48):
    return [torus_point(theta, phi, r0+(r1-r0)*i/count) for i in range(count+1)]


def build_fuel_injection_module(cols, mats, port):
    """Precámara conceptual con evacuación previa, alimentación de gas y pistón/resorte."""
    theta, phi = port["theta"], port["phi"]
    # Boquilla corta hasta el borde de plasma; el tramo a través de pared coincide con la ranura.
    nozzle_pts = _tube_points_at_minor(theta, phi, FUEL_NOZZLE_TIP_MINOR, 0.475, 72)
    nozzle = create_curve("Boquilla_Inyector_Replegada_Fuera_Del_Plasma", nozzle_pts, cols["fuel"], mats["fuel"], FUEL_NOZZLE_TUBE_RADIUS,
                          cyclic=False, resolution=1, route_kind="fuel")
    nozzle["Funcion"] = "Boquilla retraida; termina fuera de la envolvente geometrica conservadora del plasma"
    nozzle["Holgura_radial_minima_m"] = FUEL_NOZZLE_TIP_MINOR - FUEL_NOZZLE_TUBE_RADIUS - PLASMA_MINOR_MAX_BOUND
    # Cuerpo hueco de la precámara exterior.
    shell = create_oriented_hollow_tube("Camara_Esclusa_Combustible", theta, phi,
                                        FUEL_CHAMBER_START_MINOR, FUEL_CHAMBER_END_MINOR,
                                        FUEL_CHAMBER_OUTER_RADIUS, FUEL_CHAMBER_INNER_RADIUS,
                                        cols["fuel"], mats["fuel_shell"], segments=64)
    shell["Secuencia_conceptual"] = "Evacuar -> cargar gas -> aislar -> liberar pistón con resorte"
    create_annular_port_flange("Brida_Entrada_Precamara", theta, phi, 0.456,
                               outer_radius=0.031, inner_radius=0.0055, depth=0.006,
                               collection=cols["fuel"], mat=mats["flange"], segments=48)
    create_annular_port_flange("Brida_Cierre_Trasera_Precamara", theta, phi, 0.562,
                               outer_radius=0.027, inner_radius=0.010, depth=0.006,
                               collection=cols["fuel"], mat=mats["flange"], segments=48)
    # Pistón móvil y vástago; la malla representa una posición de reposo, no una simulación mecánica.
    create_oriented_cylinder("Piston_Mecanico_Inyector", theta, phi, FUEL_PISTON_MINOR,
                             radius=0.0185, depth=0.0055, collection=cols["fuel"],
                             mat=mats["fuel_piston"], segments=48)
    create_curve("Vastago_Piston_Inyector",
                 _tube_points_at_minor(theta, phi, FUEL_PISTON_MINOR+0.004, FUEL_CHAMBER_END_MINOR+0.010, 32),
                 cols["fuel"], mats["fuel_piston"], 0.0032, cyclic=False, resolution=1)
    # Resorte helicoidal en torno al eje normal (sin especificar constante elástica ni fuerza).
    e1, e2, _normal = torus_basis(theta, phi)
    spring_pts=[]
    turns=9
    steps=turns*32
    for i in range(steps+1):
        t=i/steps
        a=2.0*math.pi*turns*t
        minor=FUEL_SPRING_START_MINOR+(FUEL_SPRING_END_MINOR-FUEL_SPRING_START_MINOR)*t
        c=torus_point(theta,phi,minor)
        spring_pts.append(c + e1*(0.013*math.cos(a)) + e2*(0.013*math.sin(a)))
    spring=create_curve("Resorte_Compresion_Inyector",spring_pts,cols["fuel"],mats["fuel_spring"],0.0015,
                        cyclic=False,resolution=1)
    spring["Funcion"]="Representacion cinemática; fuerza, fatiga y energía no dimensionadas"
    # Volumen luminoso semitransparente que identifica la precámara cargada; no es una mezcla física calculada.
    gas_obj=create_uv_sphere("Volumen_Conceptual_Gas_Precamara",
                             torus_point(theta,phi,0.480),0.013,cols["fuel"],mats["fuel_gas"],segments=24,rings=12)
    gas_obj.scale=(1.0,1.0,1.8)
    gas_obj["Nota"]="Marcador visual del gas, no malla de densidad"

    # Línea de llenado con válvula y unión explícita a la botella de suministro.
    branch=FUEL_BRANCH_PHI_OFFSET
    gas_phi0=phi+0.022
    gas_phi1=phi+branch
    bottle_phi=phi+branch+0.035
    gas_points=[torus_point(theta,gas_phi0+(gas_phi1-gas_phi0)*i/64,0.515) for i in range(65)]
    create_curve("Linea_Llenado_Combustible_Precamara",gas_points,cols["fuel"],mats["fuel"],0.0022,
                 cyclic=False,resolution=1,route_kind="fuel")
    create_oriented_cylinder("Valvula_Llenado_Combustible",theta,gas_phi1,0.515,
                             radius=0.0065,depth=0.014,collection=cols["fuel"],mat=mats["fuel_shell"],segments=24)
    # El conector sale del cuerpo de válvula y llega al extremo inferior de la botella.
    bottle_connector=[torus_point(theta,gas_phi1,0.518),
                      torus_point(theta,gas_phi1,0.5185),
                      torus_point(theta,bottle_phi,0.5185)]
    create_curve("Conexion_Valvula_Llenado_A_Botella",bottle_connector,cols["fuel"],mats["fuel"],0.0022,
                 cyclic=False,resolution=1,route_kind="fuel")
    create_oriented_cylinder("Botella_Suministro_Combustible_Conceptual",theta,bottle_phi,0.548,
                             radius=0.011,depth=0.060,collection=cols["fuel"],mat=mats["fuel_shell"],segments=32)
    create_surface_text("Texto_Suministro_Combustible","SUMINISTRO GAS",theta,bottle_phi,0.581,
                        cols["labels"],mats["label_text"],size=0.008)

    # Evacuación de la precámara: la línea asciende y termina conectada al colector de escape común.
    vac_phi0=phi-0.022
    vac_phi1=phi-branch
    vac_points=[torus_point(theta,vac_phi0+(vac_phi1-vac_phi0)*i/64,0.515) for i in range(65)]
    create_curve("Linea_Evacuacion_Precamara",vac_points,cols["fuel"],mats["pump"],0.0022,
                 cyclic=False,resolution=1,route_kind="vacuum")
    create_oriented_cylinder("Valvula_Evacuacion_Precamara",theta,vac_phi1,0.515,
                             radius=0.0065,depth=0.014,collection=cols["fuel"],mat=mats["flange"],segments=24)
    lift_points=_tube_points_at_minor(theta,vac_phi1,0.515,VACUUM_HEADER_MINOR_RADIUS,64)
    create_curve("Conexion_Evacuacion_Precamara_A_Colector_Vacio",lift_points,cols["fuel"],mats["pump"],
                 0.0022,cyclic=False,resolution=1,route_kind="vacuum")
    create_surface_text("Texto_Evacuacion_Precamara","EVACUACION / VALVULA",theta,phi-branch,0.531,
                        cols["labels"],mats["label_text"],size=0.0075)
    # Elemento de disparo mecánico/seguro retenedor del pistón.
    create_oriented_box("Cerrojo_Resorte_Piston",theta,phi+0.020,0.558,
                        (0.020,0.014,0.008),cols["fuel"],mats["maintenance_latch"])
    create_surface_text("Texto_Inyector_Precamara","INYECTOR / PRECARGA + PISTON",theta,phi+0.16,0.574,
                        cols["labels"],mats["label_text"],size=0.009)
    for ob in bpy.data.objects:
        if ob.name.startswith(PREFIX) and ("Precámara" in ob.name or "Precamara" in ob.name or "Inyector" in ob.name):
            ob["Sistema"] = "Inyector de prueba conceptual; no usar como dimensionado de presión/combustible"


def build_hv_test_ignition(cols, mats, ports):
    """Dos electrodos helicoidales opuestos con pasamuros aislados en ranuras distintas."""
    sample_count=1080  # múltiplo de 6: contiene exactamente los pasos toroidales de 60° y 300°.
    # La fase se elige para que cada hélice llegue a su propio puerto en theta=pi.
    for name, spec in (("Anodo",HV_PORTS["HV_Anodo"]),("Catodo",HV_PORTS["HV_Catodo"])):
        pts=[]
        for i in range(sample_count):
            phi=2.0*math.pi*i/sample_count
            theta=HV_ELECTRODE_PHASE + spec["hand"]*POLOIDAL_TURNS_PER_TOROIDAL_TURN*phi
            pts.append(torus_point(theta,phi,spec["minor"]))
        mat=mats["hv_anode"] if name=="Anodo" else mats["hv_cathode"]
        electrode=create_curve(f"Electrodo_Helicoidal_{name}_Ensayo",pts,cols["ignition"],mat,
                               HV_ELECTRODE_RADIUS,cyclic=True,resolution=1,route_kind="hv")
        electrode["Modo"]="Ensayo de ruptura dieléctrica/ionización; no operación continua de fusión"
        electrode["Sentido_helicoidal"]="+4 poloidal/toroidal" if spec["hand"]>0 else "-4 poloidal/toroidal"
        electrode["Radio_menor_m"] = spec["minor"]
        electrode["Alimentacion"] = "Pasamuros aislado propio; no se han fijado voltaje, corriente ni energía"
        # Conductor de salida radial desde el electrodo hasta la ranura de la vasija.
        lead_pts=_tube_points_at_minor(spec["theta"],spec["phi"],spec["minor"],HV_TERMINAL_MINOR,80)
        lead=create_curve(f"Cable_HV_{name}_A_Ranura_Independiente",lead_pts,cols["ignition"],mat,
                          0.0022,cyclic=False,resolution=1,route_kind="hv")
        lead["Slot"]="HV_Anodo" if name=="Anodo" else "HV_Catodo"
        # Aisladores anulares en el paso, con agujero central para el conductor.
        for k, minor in enumerate((0.441,0.451,0.461,0.471),1):
            create_annular_port_flange(f"Aislador_Ceramico_HV_{name}_{k}",spec["theta"],spec["phi"],minor,
                                       outer_radius=0.0075,inner_radius=0.0028,depth=0.004,
                                       collection=cols["ignition"],mat=mats["hv_insulator"],segments=32)
        create_annular_port_flange(f"Brida_Terminal_HV_{name}",spec["theta"],spec["phi"],HV_TERMINAL_MINOR-0.008,
                                   outer_radius=0.013,inner_radius=0.0035,depth=0.006,
                                   collection=cols["ignition"],mat=mats["flange"],segments=32)
        create_oriented_box(f"Bloque_Terminal_HV_{name}",spec["theta"],spec["phi"],HV_TERMINAL_MINOR,
                            (0.052,0.030,0.009),cols["ignition"],mats["hv_insulator"])
        label=f"HV+ / ANODO" if name=="Anodo" else "HV- / CATODO"
        create_surface_text(f"Etiqueta_Terminal_HV_{name}",label,spec["theta"],spec["phi"],
                            HV_TERMINAL_MINOR+0.006,cols["labels"],mats["label_text"],size=0.009)
    create_surface_text("Etiqueta_Sistema_Ignicion_HV","IGNICION HV / SOLO ENSAYO",math.pi,
                        0.5*math.pi,0.490,cols["labels"],mats["label_text"],size=0.009)
    print("[V14] Ignición HV: dos electrodos helicoidales opuestos; ranuras y bornes aislados independientes.")


# =============================================================================
# 12C. ESCOTILLA, MANIPULACIÓN Y SUSTITUCIÓN DE CASSETTE
# =============================================================================
def create_curved_torus_panel(name, theta_center, phi_center, half_theta, half_phi,
                              minor_inner, minor_outer, collection, mat,
                              n_theta=28, n_phi=24):
    """Panel rectangular curvado y cerrado, definido por coordenadas toroidales."""
    verts=[]
    for minor in (minor_inner,minor_outer):
        for i in range(n_phi+1):
            phi=phi_center-half_phi+2.0*half_phi*i/n_phi
            for j in range(n_theta+1):
                theta=theta_center-half_theta+2.0*half_theta*j/n_theta
                verts.append(torus_point(theta,phi,minor))
    layer=(n_phi+1)*(n_theta+1)
    row=n_theta+1
    def ix(k,i,j): return k*layer+i*row+j
    faces=[]
    for i in range(n_phi):
        for j in range(n_theta):
            faces.append((ix(1,i,j),ix(1,i+1,j),ix(1,i+1,j+1),ix(1,i,j+1)))
            faces.append((ix(0,i,j+1),ix(0,i+1,j+1),ix(0,i+1,j),ix(0,i,j)))
    for i in range(n_phi):
        faces.append((ix(0,i,0),ix(1,i,0),ix(1,i+1,0),ix(0,i+1,0)))
        faces.append((ix(0,i+1,n_theta),ix(1,i+1,n_theta),ix(1,i,n_theta),ix(0,i,n_theta)))
    for j in range(n_theta):
        faces.append((ix(0,0,j+1),ix(1,0,j+1),ix(1,0,j),ix(0,0,j)))
        faces.append((ix(0,n_phi,j),ix(1,n_phi,j),ix(1,n_phi,j+1),ix(0,n_phi,j+1)))
    return link_mesh_object(name,verts,faces,collection,mat,smooth=False)


def create_maintenance_access(cols, mats, port):
    theta,phi=port["theta"],port["phi"]
    # Marco de junta rectangular curvada y tapa que rellena el espesor de la vasija.
    hth=port["half_theta"]+math.radians(1.2)
    hph=port["half_phi"]+math.radians(0.8)
    rframe=R_MINOR_VESSEL_OUTER+0.0007
    edges=[]
    edges.append([(theta-hth,phi-hph),(theta-hth,phi+hph)])
    edges.append([(theta+hth,phi-hph),(theta+hth,phi+hph)])
    edges.append([(theta-hth,phi-hph),(theta+hth,phi-hph)])
    edges.append([(theta-hth,phi+hph),(theta+hth,phi+hph)])
    for i,edge in enumerate(edges,1):
        pts=[]
        (t0,p0),(t1,p1)=edge
        for j in range(65):
            f=j/64.0
            pts.append(torus_point(t0+(t1-t0)*f,p0+(p1-p0)*f,rframe))
        create_curve(f"Marco_Escotilla_Mantenimiento_{i}",pts,cols["maintenance"],mats["maintenance"],
                     MAINTENANCE_FRAME_RADIUS,cyclic=False,resolution=1)
    if not MAINTENANCE_HATCH_OPEN:
        lid=create_curved_torus_panel("Tapa_Hermetica_Escotilla_Cassette_T01_C02",theta,phi,hth,hph,
                                      R_MINOR_CHAMBER-0.0010,R_MINOR_VESSEL_OUTER+0.0002,
                                      cols["maintenance"],mats["maintenance"],n_theta=32,n_phi=32)
        lid["Cierra_paso_vacio"] = True
        lid["Desmontable"] = True
        lid["Espesor_nominal_m"] = R_MINOR_VESSEL_OUTER-R_MINOR_CHAMBER+0.0012
        # Pernos grandes fuera de las zonas por donde pasan los dos puentes de bobina.
        for edge_theta in (theta-hth+math.radians(3.0), theta+hth-math.radians(3.0)):
            for edge_phi in (phi-hph+math.radians(2.0), phi, phi+hph-math.radians(2.0)):
                create_uv_sphere(f"Perno_Escotilla_{round(math.degrees(edge_theta))}_{round(math.degrees(edge_phi))}",
                                 torus_point(edge_theta,edge_phi,R_MINOR_VESSEL_OUTER+0.0012),
                                 MAINTENANCE_BOLT_RADIUS,cols["maintenance"],mats["maintenance_latch"],segments=12,rings=8)
        # Bisagras laterales compactas (el cierre se modela como brida, no como mecanismo certificado).
        for j, dth in enumerate((-math.radians(11),0.0,math.radians(11)),1):
            create_oriented_cylinder(f"Bisagra_Escotilla_{j}",theta+dth,phi+hph-0.002,
                                     R_MINOR_VESSEL_OUTER+0.0025,0.0045,0.018,
                                     cols["maintenance"],mats["maintenance_latch"],segments=20)
    # Raíles del carro de mantenimiento quedan por fuera de la envolvente de bobinas.
    rail_r=0.535
    for side, dth in (("A",-math.radians(12)),("B",math.radians(12))):
        pts=[]
        for i in range(129):
            ph=phi-math.radians(18)+math.radians(36)*i/128
            pts.append(torus_point(theta+dth,ph,rail_r))
        create_curve(f"Rail_Mantenimiento_{side}",pts,cols["maintenance"],mats["maintenance"],0.0045,
                     cyclic=False,resolution=1)
    create_oriented_box("Carro_Puente_Mantenimiento",theta,phi,0.548,(0.145,0.095,0.012),
                        cols["maintenance"],mats["maintenance_latch"])
    create_surface_text("Texto_Escotilla_Mantenimiento","ACCESO / CAMBIO DE CASSETTE",theta,
                        phi+math.radians(18),0.546,cols["labels"],mats["label_text"],size=0.009)
    if MAINTENANCE_HATCH_OPEN:
        # Ocultar el módulo instalado para mostrarlo en extracción; conservar las líneas madre.
        target_prefixes=(PREFIX+"T01_C02_",)
        for obj in bpy.data.objects:
            if obj.name.startswith(target_prefixes):
                obj.hide_set(True)
                obj.hide_render=True
        # Cable de izado y cassette de demostración, suspendidos en el exterior de la escotilla.
        cable_pts=[torus_point(theta,phi,0.548-(0.548-0.475)*i/64) for i in range(65)]
        create_curve("Cable_Izado_Cassette_Extraido",cable_pts,cols["maintenance"],mats["maintenance"],
                     0.0025,cyclic=False,resolution=1)
        sample=create_oriented_box("Cassette_Extraido_Demostracion_T01_C02",theta,phi,0.468,
                                   (0.310,0.245,0.020),cols["maintenance"],mats["backing"])
        sample["Nota"]="Pieza de demostración para modo de mantenimiento; no sustituye la cinemática de retirada"
        # Ocho placas oscuras y teselas sobre el módulo extraído.
        for row in range(2):
            for col in range(4):
                ph=phi+(col-1.5)*0.060
                th=theta+(row-0.5)*0.090
                create_oriented_box(f"Tesela_Cassette_Extraido_R{row+1}_C{col+1}",th,ph,0.482,
                                    (0.052,0.082,0.006),cols["maintenance"],mats["first_wall"])
    print("[V14] Acceso de mantenimiento: escotilla {} y cassette T01_C02.".format("ABIERTA (demo)" if MAINTENANCE_HATCH_OPEN else "cerrada"))


# =============================================================================
# 12. INFORME DE HOLGURAS Y VALIDACIONES DE CONFIGURACIÓN
# =============================================================================
def print_geometry_report(ports):
    target_outer = PLASMA_MINOR_MAX_BOUND + TILE_OUTER_OFFSET
    backing_outer = PLASMA_MINOR_MAX_BOUND + BACKING_OUTER_OFFSET
    channel_outer = PLASMA_MINOR_MAX_BOUND + CHANNEL_CENTER_OFFSET + CHANNEL_RADIUS
    header_outer = max(COOLANT_COLLECTOR_RADII.values()) + HEADER_RADIUS
    branch_outer = max(COOLANT_COLLECTOR_RADII.values()) + BRANCH_RADIUS
    pump_outer = max(
        PLASMA_MINOR_MAX_BOUND + PLENUM_FLOOR_OUTER_OFFSET + PUMP_DUCT_MOUTH_RADIUS,
        PUMP_SERVICE_MINOR_RADIUS + PUMP_DUCT_RADIUS,
    )
    # Envolvente cilíndrica de los manguitos que sobresalen de la vasija.
    ext_min, ext_max = float("inf"), float("-inf")
    for port in ports:
        if port["kind"] == "window":
            continue
        minor = port["minor_stub"]
        tube_r = HEADER_RADIUS if port["kind"] == "coolant" else PUMP_DUCT_RADIUS
        rho = R_MAJOR + minor * math.cos(port["theta"])
        ext_min = min(ext_min, rho - tube_r)
        ext_max = max(ext_max, rho + tube_r)
    report = [
        "\n========== INFORME GEOMÉTRICO REACTOR V14 ==========",
        f"Radio mayor: {R_MAJOR:.3f} m",
        f"Envolvente radial nominal del TOROIDE (solo vasija): {R_INNER_AXIS:.3f}–{R_OUTER_AXIS:.3f} m",
        f"Vasija: r_minor exterior {R_MINOR_VESSEL_OUTER:.3f} m; cámara {R_MINOR_CHAMBER:.3f} m; pared {VESSEL_WALL_THICKNESS*1000:.0f} mm",
        f"Envolvente cilíndrica de la vasija: {R_MAJOR-R_MINOR_VESSEL_OUTER:.3f}-{R_MAJOR+R_MINOR_VESSEL_OUTER:.3f} m",
        f"Plasma: radio menor máximo conservador {PLASMA_MINOR_MAX_BOUND:.3f} m",
        f"Cara posterior máxima de teselas: {target_outer:.3f} m; holgura a pared {R_MINOR_CHAMBER-target_outer:.3f} m",
        f"Respaldo máximo: {backing_outer:.3f} m; holgura a pared {R_MINOR_CHAMBER-backing_outer:.3f} m",
        f"Envolvente de canales de agua: {channel_outer:.3f} m; holgura a pared {R_MINOR_CHAMBER-channel_outer:.3f} m",
        f"Radios de agua: serpentín {CHANNEL_RADIUS*1000:.1f} mm; colector {HEADER_RADIUS*1000:.1f} mm; ramal {BRANCH_RADIUS*1000:.1f} mm",
        f"Radio de los pasajes internos de helio de bobinas: {COIL_HELIUM_CHANNEL_RADIUS*1000:.1f} mm",
        f"Envolvente de colectores: {header_outer:.3f} m; holgura a pared {R_MINOR_CHAMBER-header_outer:.3f} m",
        f"Envolvente de ramales: {branch_outer:.3f} m; holgura a pared {R_MINOR_CHAMBER-branch_outer:.3f} m",
        f"Envolvente máxima conservadora de la boca de bombeo: {pump_outer:.3f} m; holgura a pared {R_MINOR_CHAMBER-pump_outer:.3f} m",
        f"Envolvente radial de manguitos exteriores/puertos: {ext_min:.3f}–{ext_max:.3f} m (servicios auxiliares; no redefine el límite del toroide)",
        f"Bobinas: {N_HELICAL_COILS} helicoidales + {N_INCLINED_COILS} inclinadas ({N_X_PAIRS} pares); cassettes: {N_CASSETTES}; blancos: {N_CASSETTES*TILES_TOROIDAL*TILE_ROWS_PER_WING*2}",
        f"Colectores de carril fijo: {N_TRACKS*4}; puertos agua/vacio: {sum(p['kind'] != 'window' for p in ports)}; ventanas: {sum(p['kind'] == 'window' for p in ports)}",
        f"Soportes: {N_HELICAL_COILS*CASSETTES_PER_TRACK-(1 if MAINTENANCE_HATCH_OPEN else 0)} sillas H + 4 puntales H/bobina + 4 puntales por bobina X + 2 bastidores toroidales",
        f"Paquetes superconductores exteriores: {2*COIL_CONDUCTOR_RADIUS*1000:.0f}x{2*COIL_TRANSVERSE_HALF_WIDTH*1000:.0f} mm (carcasa)",
        f"Capas de bobina: núcleo {2*COIL_CORE_RADIAL_HALF_WIDTH*1000:.0f}x{2*COIL_CORE_TRANSVERSE_HALF_WIDTH*1000:.0f} mm; aislamiento {2*COIL_INSULATION_RADIAL_HALF_WIDTH*1000:.0f}x{2*COIL_INSULATION_TRANSVERSE_HALF_WIDTH*1000:.0f} mm",
        f"Bobinas helicoidales r_minor={COIL_CENTERLINE_MINOR:.3f} m; bobinas X r_minor={COIL_X_CENTERLINE_MINOR:.3f}±{COIL_X_SEPARATION_R:.3f} m",
        f"Bastidor intermedio r_minor={COIL_INNER_FRAME_MINOR:.3f}±{COIL_INNER_FRAME_HALF_THICKNESS:.3f} m; bastidor exterior={COIL_OUTER_FRAME_MINOR:.3f}±{COIL_OUTER_FRAME_HALF_THICKNESS:.3f} m",
        f"Envolvente radial calculada de bobinas X: {R_MAJOR-math.sqrt((COIL_X_CENTERLINE_MINOR+COIL_X_SEPARATION_R+COIL_CONDUCTOR_RADIUS)**2+COIL_TRANSVERSE_HALF_WIDTH**2):.3f}–{R_MAJOR+math.sqrt((COIL_X_CENTERLINE_MINOR+COIL_X_SEPARATION_R+COIL_CONDUCTOR_RADIUS)**2+COIL_TRANSVERSE_HALF_WIDTH**2):.3f} m (sin límite 0.50–1.50 m)",
        f"Conductores creados: {N_HELICAL_COILS} helicoidales + {N_INCLINED_COILS} bobinas X",
        f"Modo HV de ensayo: {HV_TEST_MODE}",
        f"CUTAWAY_VIEW={CUTAWAY_VIEW} (False = pared completa; True = corte de inspección)",
    ]
    for line in report:
        print(line)
    min_clearance = min(R_MINOR_CHAMBER - channel_outer, R_MINOR_CHAMBER - header_outer,
                        R_MINOR_CHAMBER - branch_outer, R_MINOR_CHAMBER - pump_outer)
    if min_clearance < 0.008:
        print("[AVISO] La menor holgura al volumen de cámara es menor que 8 mm; revisar parámetros.")
    coil_minor_extent = math.sqrt((COIL_X_CENTERLINE_MINOR + COIL_X_SEPARATION_R + COIL_CONDUCTOR_RADIUS)**2 + COIL_TRANSVERSE_HALF_WIDTH**2)
    coil_rmin = R_MAJOR - coil_minor_extent
    coil_rmax = R_MAJOR + coil_minor_extent
    print(f"[INFO] Envolvente de bobinas X {coil_rmin:.3f}–{coil_rmax:.3f} m; puede sobresalir porque 0.50–1.50 m solo limita el toroide.")
    if max(PLASMA_MINOR_MAX_BOUND + PLENUM_FLOOR_OUTER_OFFSET + PUMP_DUCT_MOUTH_RADIUS,
           PUMP_SERVICE_MINOR_RADIUS + PUMP_DUCT_RADIUS) >= R_MINOR_CHAMBER:
        print("[ERROR] El conducto de bombeo invade la pared. Reduce plasma, offset o radio de tubo.")
    print("[NOTA] Las holguras son cotas por envolvente. No sustituyen detección de colisiones ni validación termo-hidráulica.")
    print(f"[HOLGURA] Conducto de bombeo a la pared: {(R_MINOR_CHAMBER - PUMP_SERVICE_MINOR_RADIUS - PUMP_DUCT_RADIUS)*1000:.2f} mm.")
    print("[LEYENDA] Agua ida=azul; agua retorno=rojo; vacio/bombeo=verde; soportes helicoidales=azul; soportes X=dorado.")
    print("==================================================\n")



def validate_geometry_parameters():
    """Validaciones dimensionales previas; no sustituyen evaluación de Blender/MHD."""
    if abs(R_MINOR_CHAMBER - (R_MINOR_VESSEL_OUTER - VESSEL_WALL_THICKNESS)) > 1e-9:
        raise ValueError("Incoherencia en el espesor de la vasija.")
    if not (0.0 < FUEL_NOZZLE_TUBE_RADIUS < FUEL_NOZZLE_TIP_MINOR):
        raise ValueError("Radio/posición inválidos de la boquilla de combustible.")
    nozzle_clearance = FUEL_NOZZLE_TIP_MINOR - FUEL_NOZZLE_TUBE_RADIUS - PLASMA_MINOR_MAX_BOUND
    if nozzle_clearance < FUEL_NOZZLE_MIN_CLEARANCE:
        raise ValueError(
            f"La boquilla invade o se aproxima demasiado al plasma: holgura={nozzle_clearance*1000:.2f} mm"
        )
    backing_limit = TARGET_HALF_ANGLE + BACKING_MARGIN_ANGLE
    if math.degrees(GAP_PANEL_THETA_LO) <= math.degrees(backing_limit):
        raise ValueError("La primera pared invade el margen angular del respaldo del divertor.")
    if math.degrees(GAP_PANEL_THETA_HI) >= 90.0 - math.degrees(backing_limit):
        raise ValueError("La primera pared invade el margen angular de la pista siguiente.")
    if not (0.0 < COIL_CONDUCTOR_RADIUS and 0.0 < COIL_TRANSVERSE_HALF_WIDTH):
        raise ValueError("Sección de bobina no válida.")
    # El límite 0.50–1.50 m es solo del toroide/vasija: no se aplica a bobinas.
    helical_gap = (COIL_CENTERLINE_MINOR - COIL_CONDUCTOR_RADIUS - R_MINOR_VESSEL_OUTER)
    required_support_gap = COIL_SUPPORT_FOOT_DEPTH + COIL_HELICAL_SADDLE_SIZE[2] + 0.002
    if helical_gap < required_support_gap - 1e-9:
        raise ValueError(f"No cabe el soporte helicoidal: hueco disponible={helical_gap*1000:.2f} mm; requerido={required_support_gap*1000:.2f} mm")
    x_b_inner = COIL_X_CENTERLINE_MINOR - COIL_X_SEPARATION_R - COIL_CONDUCTOR_RADIUS
    h_outer = COIL_CENTERLINE_MINOR + COIL_CONDUCTOR_RADIUS
    interlayer_gap = x_b_inner - h_outer
    pair_radial_gap = 2.0 * COIL_X_SEPARATION_R - 2.0 * COIL_CONDUCTOR_RADIUS
    if interlayer_gap < COIL_MIN_SURFACE_CLEARANCE or pair_radial_gap < COIL_MIN_SURFACE_CLEARANCE:
        raise ValueError(f"Separación radial nominal insuficiente entre capas X/helicoidales: {interlayer_gap*1000:.2f}/{pair_radial_gap*1000:.2f} mm")
    pump_wall_clearance = R_MINOR_CHAMBER - PUMP_SERVICE_MINOR_RADIUS - PUMP_DUCT_RADIUS
    if pump_wall_clearance < 0.008:
        raise ValueError(
            f"Holgura insuficiente del conducto de bombeo a la cámara: {pump_wall_clearance*1000:.2f} mm"
        )
    if max(PLASMA_MINOR_MAX_BOUND + PLENUM_FLOOR_OUTER_OFFSET + PUMP_DUCT_MOUTH_RADIUS,
           PUMP_SERVICE_MINOR_RADIUS + PUMP_DUCT_RADIUS) >= R_MINOR_CHAMBER:
        raise ValueError("El conducto de bombeo invade la pared.")
    if not collector_service_corridors_are_clear():
        raise ValueError("Los colectores se solapan angularmente con las bandas de paneles.")
    backing_limit = TARGET_HALF_ANGLE + BACKING_MARGIN_ANGLE
    opposite_limit = 0.5 * math.pi - backing_limit
    for lane in COLLECTOR_SERVICE_OFFSETS:
        if lane - COLLECTOR_SERVICE_HALF_ANGLE <= backing_limit or lane + COLLECTOR_SERVICE_HALF_ANGLE >= opposite_limit:
            raise ValueError("Un corredor de colector invade el margen angular del cassette de divertor.")
    if CASSETTE_SPAN_FRACTION >= 1.0 or CASSETTE_SPAN_FRACTION <= 0.75:
        raise ValueError("CASSETTE_SPAN_FRACTION debe estar entre 0.75 y 1.0.")
    if MAINTENANCE_HATCH_HALF_THETA <= TARGET_HALF_ANGLE + BACKING_MARGIN_ANGLE:
        raise ValueError("La escotilla debe cubrir la sección poloidal completa del cassette.")
    return {
        "nozzle_clearance_m": nozzle_clearance,
        "helical_support_gap_m": helical_gap,
        "coil_helical_to_x_radial_gap_m": interlayer_gap,
        "coil_x_pair_radial_gap_m": pair_radial_gap,
        "pump_duct_wall_clearance_m": pump_wall_clearance,
        "first_wall_band_degrees": (math.degrees(GAP_PANEL_THETA_LO), math.degrees(GAP_PANEL_THETA_HI)),
        "panel_bands_degrees": [(math.degrees(a), math.degrees(b)) for a, b in panel_bands_avoiding_collectors()],
        "collector_lanes_degrees": [math.degrees(v) for v in COLLECTOR_SERVICE_OFFSETS],
    }



def _world_aabb(obj):
    corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    return ((min(v.x for v in corners), min(v.y for v in corners), min(v.z for v in corners)),
            (max(v.x for v in corners), max(v.y for v in corners), max(v.z for v in corners)))


def _aabb_overlaps(a, b, tol=1e-6):
    alo, ahi = a
    blo, bhi = b
    return not any(ahi[k] < blo[k] - tol or bhi[k] < alo[k] - tol for k in range(3))


def _same_support_assembly(a, b):
    """Exime solo zapata/poste/silla con el mismo identificador de soporte."""
    suffixes = ("_Zapata", "_Poste", "_Silla")
    base_a = next((a[:-len(suffix)] for suffix in suffixes if a.endswith(suffix)), None)
    base_b = next((b[:-len(suffix)] for suffix in suffixes if b.endswith(suffix)), None)
    return base_a is not None and base_a == base_b


def _strip_generated_prefix(name):
    return name[len(PREFIX):] if name.startswith(PREFIX) else name


def _support_target_coil_id(name):
    """Identifica la bobina propia a la que pertenece una silla/soporte."""
    raw = _strip_generated_prefix(name)
    if "_Silla_" not in raw:
        return None
    token = raw.split("_Silla_", 1)[0]
    if len(token) == 3 and token.startswith("H") and token[1:].isdigit():
        return ("H", token[1:])
    if len(token) == 5 and token.startswith("X") and token[1:3].isdigit() and token[3] == "_" and token[4] in "AB":
        return ("X", token[1:3], token[4])
    return None


def _is_matching_coil_pack(support_name, other_name):
    """True solo para el paquete de conductor que esa silla debe soportar."""
    key = _support_target_coil_id(support_name)
    if key is None:
        return False
    other = _strip_generated_prefix(other_name)
    if key[0] == "H":
        coil_id = key[1]
        return (f"Bobina_Helicoidal_H{coil_id}_Paquete" in other or
                f"Puente_Desmontable_Paquete_Bobina_H{coil_id}" in other)
    pos, branch = key[1], key[2]
    return (f"Bobina_X_{pos}_{branch}_Paquete" in other or
            f"Puente_Desmontable_Paquete_Bobina_X{pos}_{branch}" in other)


def _intentional_solid_contact(name_a, name_b):
    """Contactos constructivos que no deben marcarse como colisión accidental."""
    a, b = name_a, name_b
    if _same_support_assembly(a, b):
        return True
    # Zapata atornillada/apoyada en la pared de la vasija.
    if (("_Zapata" in a and "Vasija_Toroidal" in b) or
            ("_Zapata" in b and "Vasija_Toroidal" in a)):
        return True
    # La silla toca deliberadamente el paquete conductor de su bobina propia;
    # no se excluyen postes ni zapatas contra bobinas, para no ocultar errores.
    if "_Silla_" in a and _is_matching_coil_pack(a, b):
        return True
    if "_Silla_" in b and _is_matching_coil_pack(b, a):
        return True
    # Anclaje esférico apoya en el extremo del puntal/placa que fija el cassette.
    if (("_Anclaje_" in a and "_Puntal_" in b) or
            ("_Anclaje_" in b and "_Puntal_" in a)):
        return True
    # Puentes desmontables replican el mismo tramo de conductor para la vista de mantenimiento.
    if "Puente_Desmontable" in a and "Puente_Desmontable" in b:
        if a.split("_Paquete", 1)[0] == b.split("_Paquete", 1)[0]:
            return True
    return False


def _intentional_route_solid_contact(route_name, solid_name):
    """Contactos de montaje/conexión ruta-sólido documentados de forma explícita."""
    route = _strip_generated_prefix(route_name)
    solid = _strip_generated_prefix(solid_name)
    if "_Serpentin_" in route:
        base, wing = route.split("_Serpentin_", 1)
        if solid == f"{base}_Respaldo_{wing}":
            return True  # canal brazed/embebido en el reverso de su propio respaldo
    if route.startswith("Paso_"):
        port = route[len("Paso_"):]
        if solid in {f"Brida_Agua_{port}", f"Manguito_{port}", f"Conector_{port}"}:
            return True
    if route.endswith("_Enlace_Bomba_Vacio"):
        track = route.split("_", 1)[0]
        if solid == f"{track}_Carcasa_Bomba_Turbomolecular":
            return True
    if route.endswith("_Escape_Bomba"):
        track = route.split("_", 1)[0]
        if solid == f"{track}_Motor_Bomba_Turbomolecular":
            return True
    if route in {"Linea_Llenado_Combustible_Precamara", "Conexion_Valvula_Llenado_A_Botella"}:
        return solid in {"Valvula_Llenado_Combustible", "Botella_Suministro_Combustible_Conceptual"}
    if route in {"Linea_Evacuacion_Precamara", "Conexion_Evacuacion_Precamara_A_Colector_Vacio"}:
        return solid == "Valvula_Evacuacion_Precamara"
    if route == "Boquilla_Inyector_Replegada_Fuera_Del_Plasma":
        return solid in {"Brida_Entrada_Precamara", "Collar_Precámara_Inyector_Combustible"}
    if route.startswith("Cable_HV_"):
        pole = "Anodo" if "Anodo" in route else "Catodo"
        return solid.startswith((f"Aislador_Ceramico_HV_{pole}_", f"Brida_Terminal_HV_{pole}"))
    return False


def _audit_bvh_tree(obj, depsgraph):
    """Devuelve BVH y muestras en coordenadas de mundo para cualquier geometría evaluada.

    Se convierte también la geometría MESH a coordenadas de mundo, evitando mezclar
    BVH de espacio local (mallas) con BVH de espacio mundial (curvas bevel).
    Retorna (tree, sample_points, is_closed_manifold).
    """
    evaluated = obj.evaluated_get(depsgraph)
    mesh = evaluated.to_mesh()
    if mesh is None:
        return None, [], False
    try:
        verts = [evaluated.matrix_world @ vertex.co for vertex in mesh.vertices]
        polygons = [list(poly.vertices) for poly in mesh.polygons if len(poly.vertices) >= 3]
        if not verts or not polygons:
            return None, [], False

        # Muestras deterministas repartidas por la malla; sirven para detectar
        # contención completa cuando las superficies no llegan a cruzarse.
        n = len(verts)
        sample_count = min(12, n)
        if sample_count == n:
            sample_points = verts[:]
        else:
            sample_points = [verts[round(i * (n - 1) / (sample_count - 1))]
                             for i in range(sample_count)]

        # Un sólido cerrado debe tener exactamente dos caras incidentes por arista.
        edge_use = {}
        for poly in polygons:
            for a, b in zip(poly, poly[1:] + poly[:1]):
                key = (a, b) if a < b else (b, a)
                edge_use[key] = edge_use.get(key, 0) + 1
        is_closed = bool(edge_use) and all(count == 2 for count in edge_use.values())

        tree = BVHTree.FromPolygons(verts, polygons, all_triangles=False,
                                    epsilon=SOLID_AUDIT_EPSILON)
        return tree, sample_points, is_closed
    finally:
        evaluated.to_mesh_clear()


def _point_inside_closed_bvh(point, tree, max_distance=10.0):
    """Test de paridad de rayos desde un punto; requiere una malla cerrada."""
    directions = (
        Vector((1.0, 0.371390676, 0.529103)).normalized(),
        Vector((0.217391, 1.0, 0.61317)).normalized(),
        Vector((0.41917, 0.28371, 1.0)).normalized(),
    )
    votes = 0
    for direction in directions:
        origin = point.copy()
        travelled = 0.0
        hits = 0
        for _ in range(64):
            result = tree.ray_cast(origin, direction, max_distance - travelled)
            location = result[0]
            distance = result[3]
            if location is None:
                break
            hits += 1
            advance = max(1e-5, 4.0 * SOLID_AUDIT_EPSILON)
            step = max(0.0, float(distance)) + advance
            travelled += step
            if travelled >= max_distance:
                break
            origin = location + direction * advance
        if hits % 2 == 1:
            votes += 1
    return votes >= 2


def check_solid_mesh_intersections(cols):
    """Audita mallas-sólidos y curvas funcionales convertidas a tubos evaluados.

    Broad phase por AABB y narrow phase BVH. También comprueba tubos registrados
    contra sólidos; los contactos constructivos permitidos se enumeran explícitamente.
    No sustituye un B-rep, una comprobación volumétrica completa ni FEA.
    """
    if not RUN_SOLID_MESH_AUDIT:
        print("[V14] Auditoría de mallas/rutas desactivada por configuración.")
        return []

    excluded_collections = {cols["labels"].name, cols["guides"].name, cols["plasma"].name}
    candidates = []
    for obj in bpy.context.scene.objects:
        if not obj.name.startswith(PREFIX) or obj.type != "MESH":
            continue
        if any(col.name in excluded_collections for col in obj.users_collection):
            continue
        if any(token in obj.name for token in ("Texto_", "Etiqueta_", "Marcador_", "Volumen_Conceptual_Gas")):
            continue
        try:
            candidates.append((obj, _world_aabb(obj), "mesh"))
        except Exception as exc:
            print(f"[V14][AVISO] AABB no disponible para {obj.name}: {exc}")

    # Descubrir rutas por su propiedad propia RDV14_ClearanceRoute. No depende
    # del nombre ni de la ubicación de la colección: mover un objeto no lo desregistra.
    route_curve_count = 0
    for obj in bpy.context.scene.objects:
        if obj.type != "CURVE" or not bool(obj.get("RDV14_ClearanceRoute", False)):
            continue
        route_kind = str(obj.get("RDV14_ClearanceKind", ""))
        if route_kind == "coil":
            continue  # los paquetes elípticos de bobina ya se auditan como sólidos
        if route_kind not in {"cooling", "coolant", "vacuum", "fuel", "hv", "cryo_internal"}:
            audit_msg = f"{obj.name}: marca de ruta con tipo no reconocido '{route_kind}'"
            print(f"[V14][AVISO] {audit_msg}")
            continue
        try:
            bbox = _world_aabb(obj)
            pad = float(obj.get("RDV14_ClearanceRadius", 0.0)) + PIPE_CHECK_MAX_STEP
            bbox = (tuple(v - pad for v in bbox[0]), tuple(v + pad for v in bbox[1]))
            candidates.append((obj, bbox, "route"))
            route_curve_count += 1
        except Exception as exc:
            print(f"[V14][AVISO] AABB de ruta no disponible para {obj.name}: {exc}")

    if len(candidates) > SOLID_AUDIT_MAX_OBJECTS:
        print(f"[V14][INCONCLUSO] {len(candidates)} objetos superan el límite {SOLID_AUDIT_MAX_OBJECTS}; no se declara ausencia de colisiones.")
        return []

    candidates.sort(key=lambda item: item[1][0][0])
    active = []
    possible_pairs = []
    for item in candidates:
        obj, bbox, kind = item
        active = [old for old in active if old[1][1][0] >= bbox[0][0] - 1e-6]
        for other in active:
            other_obj, other_bbox, other_kind = other
            if not _aabb_overlaps(other_bbox, bbox):
                continue
            if kind == "route" and other_kind == "route":
                continue  # ruta-ruta ya está cubierta por check_all_pipe_clearances()
            if kind == "mesh" and other_kind == "mesh":
                assembly_a = other_obj.get("RDV14_CoilAssembly")
                assembly_b = obj.get("RDV14_CoilAssembly")
                if assembly_a and assembly_a == assembly_b:
                    continue  # núcleo/aislamiento/carcasa y puntales de la misma bobina
                frame_a = other_obj.get("RDV14_SupportFrameId")
                frame_b = obj.get("RDV14_SupportFrameId")
                if frame_a and frame_a == frame_b:
                    continue  # puntal/collar unido intencionadamente a su bastidor toroidal
                if _intentional_solid_contact(other_obj.name, obj.name):
                    continue
            else:
                route_obj = obj if kind == "route" else other_obj
                solid_obj = other_obj if kind == "route" else obj
                if (route_obj.get("RDV14_ClearanceKind") == "cryo_internal"
                        and route_obj.get("RDV14_CoilAssembly")
                        and route_obj.get("RDV14_CoilAssembly") == solid_obj.get("RDV14_CoilAssembly")):
                    continue  # pasaje de helio/toma dentro de la carcasa de su bobina
                route_frame = route_obj.get("RDV14_SupportFrameId")
                solid_frame = solid_obj.get("RDV14_SupportFrameId")
                if route_frame and route_frame == solid_frame:
                    continue  # toma de helio que entra intencionadamente en el bastidor
                if _intentional_route_solid_contact(route_obj.name, solid_obj.name):
                    continue
            possible_pairs.append((other_obj, obj, other_kind, kind))
        active.append(item)

    if len(possible_pairs) > SOLID_AUDIT_MAX_PAIRS:
        print(f"[V14][INCONCLUSO] {len(possible_pairs)} pares superan el límite {SOLID_AUDIT_MAX_PAIRS}; no se declara ausencia de colisiones.")
        return []

    depsgraph = bpy.context.evaluated_depsgraph_get()
    bvh_cache = {}
    intersections = []
    containments = []
    audit_errors = []
    containment_checks = 0
    containment_check_limit = 6000
    for a, b, kind_a, kind_b in possible_pairs:
        try:
            if a.name not in bvh_cache:
                bvh_cache[a.name] = _audit_bvh_tree(a, depsgraph)
            if b.name not in bvh_cache:
                bvh_cache[b.name] = _audit_bvh_tree(b, depsgraph)
            tree_a, samples_a, closed_a = bvh_cache[a.name]
            tree_b, samples_b, closed_b = bvh_cache[b.name]
            if tree_a is None or tree_b is None:
                audit_errors.append((a.name, b.name, "BVH no disponible/geometría vacía"))
                continue
            overlaps = tree_a.overlap(tree_b)
            if overlaps:
                intersections.append((a.name, b.name, len(overlaps)))
                continue

            # BVH.overlap detecta cruces de superficies, pero no un sólido completamente
            # contenido dentro de otro. Para pares con AABB solapadas se hace una prueba
            # de paridad con muestras de vértices contra el BVH del sólido cerrado.
            contained = False
            frame_shell_pair = bool(a.get("RDV14_SupportFrameShell") or b.get("RDV14_SupportFrameShell"))
            if closed_b and not frame_shell_pair:
                for sample in samples_a[:4]:
                    containment_checks += 1
                    if containment_checks > containment_check_limit:
                        break
                    if _point_inside_closed_bvh(sample, tree_b):
                        contained = True
                        break
            if not contained and closed_a and containment_checks <= containment_check_limit and not frame_shell_pair:
                for sample in samples_b[:4]:
                    containment_checks += 1
                    if containment_checks > containment_check_limit:
                        break
                    if _point_inside_closed_bvh(sample, tree_a):
                        contained = True
                        break
            if contained:
                containments.append((a.name, b.name))
        except Exception as exc:
            audit_errors.append((a.name, b.name, str(exc)))
            print(f"[V14][AVISO] No se pudo comprobar {a.name} / {b.name}: {exc}")

    print("\n========== AUDITORÍA BVH V14 ==========")
    print(f"Mallas sólidas: {sum(1 for _,_,kind in candidates if kind == 'mesh')}")
    print(f"Curvas funcionales convertidas a tubos: {route_curve_count}")
    print(f"Pares AABB candidatos: {len(possible_pairs)}")
    print(f"Pares con intersección de superficies: {len(intersections)}")
    print(f"Pares con contención detectada: {len(containments)}")
    print(f"Pruebas de contención realizadas: {min(containment_checks, containment_check_limit)}")
    for a, b, count in intersections[:SOLID_AUDIT_MAX_REPORTED]:
        print(f"  [REVISAR SUPERFICIE] {a} <> {b}: {count} pares de caras intersectados")
    for a, b in containments[:SOLID_AUDIT_MAX_REPORTED]:
        print(f"  [REVISAR CONTENCIÓN] {a} <> {b}: una muestra queda dentro de un sólido cerrado")
    if len(intersections) > SOLID_AUDIT_MAX_REPORTED:
        print(f"  ... {len(intersections)-SOLID_AUDIT_MAX_REPORTED} pares de superficie adicionales omitidos")
    if len(containments) > SOLID_AUDIT_MAX_REPORTED:
        print(f"  ... {len(containments)-SOLID_AUDIT_MAX_REPORTED} pares de contención adicionales omitidos")
    if containment_checks > containment_check_limit:
        print("[INCONCLUSO] Se alcanzó el límite de pruebas de contención; quedan pares sin prueba de anidamiento.")
    if intersections or containments:
        print("[REVISAR] Hay cruces o contenciones que requieren inspección.")
    elif audit_errors:
        print(f"[INCONCLUSO] {len(audit_errors)} pares no pudieron comprobarse.")
    elif containment_checks > containment_check_limit:
        print("[INCONCLUSO] No se declaran todas las colisiones libres porque el límite de contención fue alcanzado.")
    else:
        print("[OK] No se detectaron cruces de superficie ni contenciones en los pares examinados.")
    print("[NOTA] El test de contención usa paridad de rayos sobre mallas cerradas; no sustituye B-rep, tolerancias, FEA ni inspección visual en Blender.")
    print("=======================================\n")
    return intersections




def export_v15_magnetic_geometry(filepath=None, sample_count=V15_MAGNETIC_EXPORT_SAMPLES):
    """Exporta las 20 centrolineas exactas V14 para el motor Biot-Savart V15.

    Las corrientes no se deducen de los materiales ni de la escena CAD. Se guardan
    como ratios de referencia explícitos, separados de la geometría. El signo de
    corriente es positivo siguiendo el orden de puntos exportado.
    """
    filepath = filepath or V15_MAGNETIC_EXPORT_PATH
    output_path = bpy.path.abspath(filepath)
    specs = build_coil_centerline_samples(sample_count=int(sample_count))
    coils = []
    for name, points, _radius in specs:
        if name.startswith("H"):
            group = "H"
        else:
            position = int(name[1:3]) - 1
            layer = name[-1]
            family = "eq" if position % 2 == 0 else "pol"
            group = f"X{family}{layer}"
        coils.append({
            "name": name,
            "group": group,
            "current_ratio_Ih": float(V15_CURRENT_REFERENCE_RATIOS[group]),
            "point_order": "H: phi increasing; X: theta increasing",
            "xyz_m": [[float(p.x), float(p.y), float(p.z)] for p in points],
        })

    track_phases = [THETA_PHASE + 2.0 * math.pi * k / N_TRACKS for k in range(N_TRACKS)]
    data = {
        "schema": "reactor-v15-magnetic-geometry-1",
        "source": {
            "repository": "wili-201027/TO_HELP",
            "branch": "v15-magnetostatic-integration",
            "source_file": "reactor_pared_divertor_v14.py",
            "source_blob_sha": "9356deb9d9614cb8789b6df2c3b24640cf072a62",
            "geometry_basis": "Exact V14 parametric centerline equations; normal installed configuration closes maintenance bridges."
        },
        "model": {
            "regime": "stationary vacuum magnetostatics",
            "dB_dt": 0.0,
            "plasma_current_included": False,
            "hv_electrodes_included": False,
            "confinement_assumed": False,
            "current_sign_convention": "Positive current follows the exported ordered centerline points."
        },
        "geometry": {
            "R_major_m": R_MAJOR,
            "plasma_radius_base_m": R_PLASMA_BASE,
            "plasma_helical_amplitude_m": AMP_HELICAL,
            "plasma_quadrupolar_amplitude_m": AMP_QUADRUPOLAR,
            "helical_count": N_HELICAL_COILS,
            "helical_minor_radius_m": COIL_CENTERLINE_MINOR,
            "helical_poloidal_turns_per_toroidal_turn": POLOIDAL_TURNS_PER_TOROIDAL_TURN,
            "x_count": N_INCLINED_COILS,
            "x_positions": 8,
            "x_tilt_rad": COIL_X_TILT,
            "x_layer_A_minor_radius_m": COIL_X_CENTERLINE_MINOR + COIL_X_SEPARATION_R,
            "x_layer_B_minor_radius_m": COIL_X_CENTERLINE_MINOR - COIL_X_SEPARATION_R,
            "x_equatorial_positions_zero_based": [0, 2, 4, 6],
            "x_polar_positions_zero_based": [1, 3, 5, 7],
            "track_phases_rad": track_phases,
            "target_half_angle_rad": TARGET_HALF_ANGLE,
            "pump_slot_half_angle_rad": PUMPING_GAP_HALF,
            "cassette_count_per_track": CASSETTES_PER_TRACK,
            "cassette_span_fraction": CASSETTE_SPAN_FRACTION,
            "target_face_offset_from_plasma_m": TILE_FACE_OFFSET,
            "pump_floor_inner_offset_from_plasma_m": PLENUM_FLOOR_INNER_OFFSET,
            "first_wall_inner_offset_from_plasma_m": FIRST_WALL_INNER_OFFSET,
            "first_wall_panel_bands_rad": [list(x) for x in panel_bands_avoiding_collectors()],
            "collector_service_offsets_rad": list(COLLECTOR_SERVICE_OFFSETS),
            "collector_service_half_angle_rad": COLLECTOR_SERVICE_HALF_ANGLE,
            "vessel_inner_minor_radius_m": R_MINOR_CHAMBER
        },
        "current_reference_ratios": dict(V15_CURRENT_REFERENCE_RATIOS),
        "currents_note": "Ratios adopted from the corrected study as a reference scenario; sign convention and net X toroidal contribution must be checked by the field audit.",
        "coil_count": len(coils),
        "sample_count_per_coil": int(sample_count),
        "coils": coils,
    }
    with open(output_path, "w", encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
    print(f"[V15] Geometria magnetica exportada: {output_path}")
    print("[V15] Export incluye 4 H + 16 X; HV no forma parte del modelo magnetico.")
    return output_path

# =============================================================================
# 13. CONSTRUCCIÓN COMPLETA
# =============================================================================
def build_reactor_divertor_v15():
    validation = validate_geometry_parameters()
    print("[V14] Validaciones dimensionales: ", validation)
    check_coil_clearances()
    clean_old_system()
    cols=create_collections()
    mats=make_materials()
    ports=build_port_specs()

    # Vasija con aberturas paramétricas: agua, vacío, inyector, HV, mantenimiento y ventana.
    create_vessel_shell(cols["vessel"],mats,ports)
    for port in ports:
        if port["kind"] in ("window","maintenance"):
            continue
        create_port_collar("Collar_"+port["name"],port,cols["feedthroughs"],mats)
        if port["kind"] == "coolant":
            create_port_tube(port,cols["feedthroughs"],mats)
            create_coolant_connector_and_label(port,cols,mats)
    window_port=next(p for p in ports if p["kind"]=="window")
    create_window_assembly(window_port,cols,mats)

    legacy_plasma = bpy.data.objects.get("Plasma_Volume_Combined")
    if CREATE_PLASMA_REFERENCE and not HV_TEST_MODE:
        if legacy_plasma is None or REMOVE_LEGACY_PLASMA_OBJECT:
            create_plasma_reference(cols["plasma"],mats["plasma"])
        else:
            print("[V14] Se conserva Plasma_Volume_Combined y no se crea una segunda superficie de plasma.")

    total_cassettes=0
    for track in range(N_TRACKS):
        phase=THETA_PHASE+2.0*math.pi*track/N_TRACKS
        build_helical_headers(track,phase,cols,mats)
        build_pumping_plenum(track,phase,cols,mats)
        build_pump_duct(track,phase,cols,mats)
        create_vacuum_pump_module(track,phase,cols,mats)
        for segment in range(CASSETTES_PER_TRACK):
            phi_center,phi_min,phi_max=cassette_phi_bounds(segment)
            build_cassette(track,segment,phase,cols,mats)
            build_cassette_cooling(track,segment,phase,phi_min,phi_max,cols,mats)
            create_supports(track,segment,phase,phi_center,cols,mats)
            total_cassettes+=1

    build_first_wall_blanket_shield(cols,mats,ports)
    fuel_port=next(p for p in ports if p["kind"]=="fuel")
    build_fuel_injection_module(cols,mats,fuel_port)
    build_vacuum_exhaust_system(cols,mats)
    if HV_TEST_MODE:
        build_hv_test_ignition(cols,mats,ports)
        print("[V14] Modo HV de ensayo activo; no se crea el plasma de referencia.")
    else:
        print("[V14] Electrodos HV omitidos en modo normal; activar HV_TEST_MODE solo para inspección de ensayo.")

    build_actual_coil_conductors(cols,mats)
    build_coil_supports(cols,mats)
    maintenance_port=next(p for p in ports if p["kind"]=="maintenance")
    create_maintenance_access(cols,mats,maintenance_port)

    if EXPORT_V15_MAGNETIC_GEOMETRY:
        export_v15_magnetic_geometry()

    bpy.ops.object.select_all(action="DESELECT")
    bpy.context.view_layer.update()
    check_all_pipe_clearances()
    check_solid_mesh_intersections(cols)
    print_geometry_report(ports)
    print("[V14] Construcción completada.")
    print(f"[V14] Cassettes de divertor: {total_cassettes}; primera pared/blanket/blindaje en sectores disponibles.")
    print("[V14] Bobinas: carcasa 48x64 mm + aislamiento + haz superconductor; doble paso de helio con puente U y dos interfaces por bobina.")
    print("[V14] Soportes: bastidor intermedio para H/XB y bastidor exterior para XA; puntales cortos en huecos dedicados.")
    print("[V14] Los radios 0.50–1.50 m se aplican solo al toroide; la envolvente de bobinas se informa por separado.")
    print("[V14] Inyector: línea de combustible separada y evacuación de precámara conectada al colector común de vacío.")
    if HV_TEST_MODE:
        print("[V14] Electrodos opuestos de ensayo HV presentes; sin fuente ni parámetros eléctricos de operación.")
    print("[V14] Escotilla de servicio para T01_C02; modo abierto separa los puentes H01/X02B y muestra una pieza de demostración.")
    print("[V14] Primera pared, blanket y blindaje cubren las bandas libres entre pistas, con huecos de portales; cobertura y espesores siguen siendo conceptuales.")
    print("[V14] La boquilla queda retraída fuera de la envolvente geométrica conservadora del plasma; el pistón/resorte sigue siendo una maqueta no dimensionada.")
    print("[V14] La ruta de evacuación de la precámara termina en el colector anular de la planta de vacío.")


if __name__ == "__main__":
    build_reactor_divertor_v15()
