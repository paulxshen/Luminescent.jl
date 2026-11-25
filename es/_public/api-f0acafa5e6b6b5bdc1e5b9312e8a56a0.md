# Referencia de la API
## Flujo de trabajo
### mesh
```python
lumi.mesh(path, component, layer_stack, materials, zmin, zmax, tolerance=0.001)
```
Escribe mallas .obj del `component` de gdsfactory extruido según `layer_stack`, recortado a [zmin, zmax], en la carpeta de simulación `path`. Llámese antes de `setup`.
Los polígonos se simplifican dentro de `tolerance` [um] antes de la extrusión; 0 lo desactiva.
También guarda el componente como layout.gds en `path`.

### setup
```python
lumi.setup(path, sources, monitors, nres=6, wls=None, wl=None, f=None, fs=None, c0=None, targets=[], designs=[], approx_2D_mode=None, **kwargs)
```
Escribe el problema de simulación en la carpeta `path`. Llámese primero a `mesh` y después a `solve`.

Argumentos:

- wl y/o wls: longitud de onda central y longitudes de onda en las que calcular los resultados, o
- f y/o fs con c0: frecuencia central y frecuencias, con la velocidad de la luz c0 en sus unidades, p. ej. 300 [mm GHz]. En ese caso, los resultados se ordenan por frecuencia
- sources, monitors: de `source`, `monitor`, `gf_source`, `gf_monitor`
- frames: de `plane_frame`, `box_frame`, `spherical_frame`, `points_frame`, `cage_frame`
- materials: diccionario de nombre a `material`, p. ej. de `MATERIALS`
- length_unit: unidad de longitud en m, p. ej. 1e-6, requerida por las entradas dispersivas de `MATERIALS`
- background: nombre del material que rellena el espacio vacío
- bbox: caja de simulación [[xmin, ymin, zmin], [xmax, ymax, zmax]]. Las estructuras cortadas en su borde continúan sin cambios más allá de él (dentro de la PML)
- boundaries: por eje, "PML", "PEC", "PMC", "periodic" o un par (lado mínimo, lado máximo), p. ej. ("PEC", "PML")
- k_bloch: vector de onda de Bloch [kx, ky, kz] en rad por unidad de longitud en los ejes periódicos, p. ej. de `k_bloch`. Los campos adquieren una fase e^(i k L) a lo largo de cada periodo. Un valor no nulo hace que los campos sean complejos (2x memoria)
- nres: puntos de malla por longitud de onda en el material. 6 es un buen punto de partida para fotónica; auméntese para comprobar la convergencia
- tmax: tiempo máximo de simulación [periodos]. Necesario para resolver; si no se establece, solo se ejecuta la configuración
- shutoff: detener antes de `tmax` una vez que la energía del campo decae a esta fracción. Requiere `tmax`. Al menos 1e-3 cuando nres <= 10
- PML_decay: atenuación del campo en la PML, escalar o por eje
- courant: paso temporal relativo al máximo estable. Si la ejecución se vuelve inestable, se reintenta con subpixel_smoothing 0.75, luego 0.5, luego courant 0.75 con subpixel_smoothing 0.5 y, finalmente, 0 (los valores establecidos nunca se aumentan)
- subpixel_smoothing: promediado subpíxel de las interfaces entre materiales: 1 anisótropo (preciso), 0 promedio volumétrico escalar (estable), valores intermedios una combinación de ambos; negativo para ninguno (escalonado)
- views: de `view`, para la animación
- shot_every: guardar un fotograma de la animación cada _ periodos
- full_fields: guardar instantáneas completas 3D de todas las componentes del campo en lugar de solo cortes de los campos de las vistas, de modo que `movie(path, views=...)` pueda volver a renderizar nuevos cortes y campos
- device: "gpu" (por defecto) o "cpu"
- approx_2D_mode: aproximación 2D
- designs, targets, optimizer: diseño inverso, de `design`, `target`, `optimizer`

### solve
```python
lumi.solve(path, fields=False, gzip=False, run=True)
```
Ejecuta la simulación en la carpeta `path` preparada por `setup`. Escribe solution.json.
`run=False`: solo configura (malla, modos de los puertos) y se detiene antes del avance temporal.
Se ejecuta en el servidor Colab de Lumi (véase `remote_solve`; `fields` y `gzip` se le pasan).

### load
```python
lumi.load(path, show=True)
```
Carga el diccionario de la solución desde `path`. Los resultados de los monitores están en sol["monitors"][name], p. ej. ["modal"][número de modo]["+" o "-"] como array sobre las longitudes de onda (o frecuencias si se configuró con fs). Los arrays en el dominio de la frecuencia tienen la longitud de onda como último eje: ["fields"][componente], p. ej. [x, y, longitud de onda], ["flux"], ["power"], ["farfield"] y ["directivity"] [theta, phi, longitud de onda], y en el monitor de tipo cage ["flux"][cara]. Los modos volumétricos de los monitores de caja no tienen dirección: ["modal"][número de modo] es la amplitud modal media después de t0, y |a|^2 la energía almacenada.

### movie
```python
lumi.movie(path, views=None, fps=None, src='frames')
```
Muestra la película de la simulación en `path`, renderizada por `solve`.
`views` (sobrescribe las vistas de setup) o `fps` la vuelven a renderizar a partir de instantáneas de campo, que las resoluciones en la nube no conservan.
También se vuelve a renderizar si las instantáneas son más recientes que la película.
`src`: "frames" para instantáneas temporales, "history" para la instantánea de pico por iteración de optimización.

### view
```python
lumi.view(field, name='', prop='n', x=None, y=None, z=None, field_color_intensity=1.0, show_grid=False, show_box=True, labels='auto', mirrors=[])
```
Vista de película: cortes de `field` superpuestos sobre la propiedad del material `prop`, para `setup(views=...)`.

Args:

- field: componente de campo, p. ej. "Ey", "Hz"
- name: nombre del archivo de la película; por defecto, el campo y las posiciones de los cortes
- prop: propiedad del material subyacente, p. ej. "n", "epsilon"
- x, y, z: posición del corte a lo largo de ese eje, "half" para el centro de la caja de simulación, o None para no hacer corte. Un panel por corte
- field_color_intensity: escala la saturación del color del campo
- show_grid: dibuja las líneas de la malla
- show_box: dibuja el contorno de la caja de simulación
- labels: "auto" para etiquetar fuentes y monitores, o un dict de nombre a etiqueta
- mirrors: refleja respecto a las fronteras de simetría para la visualización, p. ej. ["z-"]

### gui
```python
lumi.gui(path=None, port=None, open=True)
```
GUI en el navegador de las ejecuciones de simulación: geometría 3D, fuentes, monitores y fronteras; parámetros impresos y registro; película de campos y pestañas de resultados (S-parameters, modos resueltos, campos, trazas temporales, datos sin procesar de los monitores) tras la resolución. El panel de código de setup se puede editar: Save escribe el code.py de la ejecución y el script que ejecutó setup (la versión anterior se conserva como .py.bak); Run vuelve a ejecutar ese script, con la salida en el panel de registro. El selector de ejecuciones lista las carpetas junto a `path` y en ./runs, ./run. Sin `path`, sigue la ejecución más reciente, también al cambiar una ventana abierta (opción conmutable en la GUI). Se actualiza a medida que cambian los archivos.
Una llamada posterior reutiliza una ventana abierta si muestra `path` o una ejecución en la misma carpeta, o si sigue la ejecución más reciente (cambiándola y añadiendo ejecuciones); en caso contrario, inicia otro servidor en el siguiente puerto. Devuelve el enlace de localhost.
`open`: la abre en el navegador, salvo que ya haya abierta una ventana reutilizada (en sesiones remotas de VS Code, también redirige el puerto). Por defecto: True.

## Materiales
### material
```python
lumi.material(n=None, epsilon=None, mu=1, sigma=0, chi2=0, chi3=0, dispersion=[], placeholder=False, dx=None)
```
Define las propiedades del material en unidades relativas

Args:

- epsilon: permitividad relativa real; el valor instantáneo (a frecuencia infinita) cuando se establece `dispersion`. Use `sigma` o `dispersion` para las pérdidas
- mu: permeabilidad relativa
- sigma: conductividad (multiplique por `377 * (characteristic wavelength in m)` para convertir de unidades métricas (`S/m, 1 /(Ω m)`) a unidades relativas)
- chi2: susceptibilidad no lineal de segundo orden de Pockels
- chi3: susceptibilidad no lineal de tercer orden de Kerr
- dispersion: lista de términos de `lorentzian`, `drude` o `gain`, p. ej. `[lumi.lorentzian(f, gamma, sigma)]`. Usamos la convención de Meep https://meep.readthedocs.io/en/latest/Materials/#material-dispersion y unidades relativas.
- n: densidad de malla real relativa al vacío; por defecto, la raíz cuadrada de la permitividad relativa. Establece epsilon si no se proporciona epsilon
- dx: espaciado máximo de la malla en este material, en unidades de longitud de la geometría. Aumenta la densidad de malla a `wavelength / (dx * nres)` en `setup`, nunca la reduce. Requiere n o epsilon

### PEC
```python
lumi.PEC(n=1, dx=None)
```
Material conductor eléctrico perfecto para modelar metales

### placeholder_material
```python
lumi.placeholder_material(n=1, dx=None)
```
El material de marcador de posición contribuye al mallado pero no a las propiedades físicas; p. ej., se usa para la región de diseño

### lorentzian
```python
lumi.lorentzian(f, gamma, sigma)
```
Término de dispersión lorentziana para `material(dispersion=[...])`: frecuencia de resonancia `f`, amortiguamiento `gamma`, intensidad `sigma` (convención de Meep, unidades relativas)

### drude
```python
lumi.drude(f, gamma, sigma)
```
Término de dispersión de Drude para `material(dispersion=[...])`: frecuencia de tipo plasma `f`, amortiguamiento `gamma`, intensidad `sigma` (convención de Meep, unidades relativas)

### gain
```python
lumi.gain(f, gamma, sigma)
```
Término de ganancia lorentziano para `material(dispersion=[...])`: como `lorentzian` con intensidad negativa

## Fuentes y monitores
### source
```python
lumi.source(name, origin, basis=None, frame=None, dim=None, mirror=None, wl=None, f=None, bandwidth=None, duration=None, modenums=[0])
```
Fuente de pulso gaussiano modulado con atenuación espectral gradual en los límites del ancho de banda. Se convierte en onda continua si `bandwidth=0`.

Argumentos:

- name: debe coincidir con un nombre de puerto, p. ej. "o1"
- origin: posición [x, y, z], o coordenada a lo largo de `dim` para marcos `face`
- basis: filas x local, y local, normal saliente. La fuente inyecta de forma bidireccional (a lo largo de ±normal)
- frame: nombre del `plane_frame` o `face` que contiene el perfil modal
- dim: normal saliente a lo largo del eje ±1, ±2, ±3 (x, y, z) en lugar de `basis`
- mirror: 'x' o 1, 'y' o 2: refleja el marco y los modos respecto a x=0 o y=0 locales, p. ej. para un puerto en el lado opuesto de un marco simétrico. None: sin reflexión
- wl o f: longitud de onda o frecuencia, especifique una de ellas (no ambas). f requiere c0 en setup
- bandwidth: en las mismas unidades que lo anterior
- modenums: lista de números de modo a excitar

### monitor
```python
lumi.monitor(name, origin, basis=None, frame=None, dim=None, mirror=None, t0=0, x=None, y=None, z=None, direction='+', layer=None)
```
Monitor que registra los resultados de `frame` en `origin` en `basis`.

Argumentos:

- name: debe coincidir con un nombre de puerto, p. ej. "o1". Único entre los monitores
- origin: posición [x, y, z], o coordenada a lo largo de `dim` para marcos `face`
- basis: filas x local, y local, normal saliente. Los resultados "+" se propagan a lo largo de la normal, los "-" en sentido contrario
- frame: nombre del marco
- dim: normal saliente a lo largo del eje ±1, ±2, ±3 (x, y, z) en lugar de `basis`
- mirror: 'x' o 1, 'y' o 2: refleja el marco y los modos respecto a x=0 o y=0 locales. None: sin reflexión
- t0: comienza a registrar señales en el dominio del tiempo en t0 [periodos]

### gf_source
```python
lumi.gf_source(port, c, layer_stack, offset=0, name=None, **kwargs)
```
`source` en el puerto `port` de gdsfactory del componente `c`, con `offset` a lo largo de la normal saliente. `name` toma por defecto el nombre del puerto. Los kwargs se pasan a `source`.

### gf_monitor
```python
lumi.gf_monitor(port, c, layer_stack, offset=0, name=None, **kwargs)
```
`monitor` en el puerto `port` de gdsfactory del componente `c`. `name` toma por defecto el nombre del puerto. Los kwargs se pasan a `monitor`.

### gaussian_mode
```python
lumi.gaussian_mode(D, n, dx, wl, theta=0, half=None, Dlim=None)
```
Modo de haz gaussiano inclinado para `plane_frame(modes=...)`, p. ej. fibra.

Argumentos:

- D: diámetro del campo modal (intensidad 1/e^2)
- n: índice de refracción del medio
- dx: muestreo
- wl: longitud de onda en el vacío
- theta: ángulo de inclinación [rad] en el medio, en el plano xz local
- half: "y+" para medio marco y>=0 (plano de simetría en y=0)
- Dlim: ancho total del marco, por defecto 1.5*D

## Marcos
### plane_frame
```python
lumi.plane_frame(name, start, stop, fs_mode_solver=[], wls_mode_solver=[], modes='auto', nmodes=1, metallic_boundaries=None, boundaries=('PML', 'PML'), time_domain=None, frequency_domain=None)
```
Marco de puerto plano entre `start` y `stop` locales (2 coordenadas en 3D). Soluciones modales en las frecuencias o longitudes de onda indicadas; las fuentes y los monitores usan (o interpolan) los modos en las longitudes de onda más cercanas. Resultados: modales, flujo, potencia. Registra resultados en el dominio de la frecuencia por defecto, o señales en el dominio del tiempo si `time_domain=True` (o `frequency_domain=False`).

Argumentos:

- fs_mode_solver o wls_mode_solver: frecuencias (requiere c0 en setup) o longitudes de onda en las que resolver los modos; especifique una de ellas (no ambas)
- modes: "auto" para resolverlos; de lo contrario, campos precalculados, p. ej. [[{'Ex': np_array, 'Hy': ...}, ...], ...]
- nmodes: número de modos a resolver
- metallic_boundaries: contornos PEC para el solucionador modal: nombres de capas o bordes del marco "x-", "x+", "y-", "y+"

### points_frame
```python
lumi.points_frame(name, points, sensor=None, time_domain=None, frequency_domain=None)
```
Puntos en la base del monitor. Registra resultados en el dominio de la frecuencia por defecto, o señales en el dominio del tiempo si `time_domain=True` (o `frequency_domain=False`).

Argumentos:

- sensor: "voltage": ∫E·dl a lo largo de la polilínea `points`. "current": ∮H·dl alrededor del lazo cerrado `points`. None: campos sin procesar.

### box_frame
```python
lumi.box_frame(name, start, stop, dr, modes=None, wavelengths=[], time_domain=None, frequency_domain=None)
```
Marco de volumen en forma de caja entre `start` y `stop` locales (3 coordenadas en 3D), muestreado cada `dr` (escalar o por eje). Resultados: campos en la base del monitor. Registra resultados en el dominio de la frecuencia por defecto, o señales en el dominio del tiempo si `time_domain=True` (o `frequency_domain=False`).

- modes: modos de volumen opcionales, p. ej., modos de resonador [[{'Ex': np_array, ..., 'Hz': ...}, ...], ...] por longitud de onda, muestreados sobre la caja (sustituye a `dr`). En ese caso, los resultados también incluyen modal: proyección media, después del t0 del monitor, sobre cada modo normalizado a energía almacenada unitaria, de modo que |a|^2 es la energía en el modo. Sin dirección: usar en los objetivos como p. ej. "o2@0"
- wavelengths: longitudes de onda de `modes`. Si está vacío, se usan los mismos modos en todas las longitudes de onda

### cage_frame
```python
lumi.cage_frame(name, start, stop, dr, theta=3.141592653589793, phi=6.283185307179586, dtheta=np.float64(0.08726646259971647), dphi=np.float64(0.08726646259971647), time_domain=None, frequency_domain=None)
```
Jaula: superficie cerrada de la caja entre `start` y `stop` locales (solo 3D), con caras muestreadas cada `dr` (escalar o por eje). Resultados: flujo por cara (hacia fuera), potencia, intensidad de campo lejano y directividad en la malla (theta, phi). theta, phi: ángulo máximo (rango desde 0) o [mín, máx]. Registra resultados en el dominio de la frecuencia por defecto, o señales en el dominio del tiempo si `time_domain=True` (o `frequency_domain=False`).

### spherical_frame
```python
lumi.spherical_frame(name, r, theta=3.141592653589793, phi=6.283185307179586, dtheta=np.float64(0.08726646259971647), dphi=np.float64(0.08726646259971647), time_domain=None, frequency_domain=None)
```
Esfera (círculo en 2D) de radio `r` alrededor del origen del monitor, en la base del monitor. theta, phi: ángulo máximo (rango desde 0) o [mín, máx]. Resultados: flujo (radial), potencia y, si es cerrada (theta y phi completos), intensidad de campo lejano y directividad en la malla (theta, phi). Registra resultados en el dominio de la frecuencia por defecto, o señales en el dominio del tiempo si `time_domain=True` (o `frequency_domain=False`).

## Diseño inverso
### design
```python
lumi.design(name, layer, lmin, fill_material, void_material, component, layer_stack, materials, symdims=[], repdims=[3], anchordims=[], dependent=False, topopt=True, contour_layer=None)
```
Región de diseño que se va a optimizar

Argumentos:

- name: nombre del diseño, p. ej., "d1"
- layer: capa de marcador de posición para la región de diseño
- lmin: tamaño mínimo de característica
- fill_material: p. ej., 'Si'
- void_material: p. ej., 'SiO2'
- symdims: p. ej., ["x","y","diagonal"]
- component, layer_stack, materials: como en setup
- repdims: p. ej., [3] para una región de diseño planar extruida a lo largo de z
- anchordims: dimensiones con signo, p. ej., [-1, 2], cuyos bordes de diseño (x-, y+) están cerrados en lugar de abiertos
- topopt: permitir cambios de topología
- contour_layer: capa opcional cuyos polígonos enmascaran la región de diseño; solo se optimizan los píxeles del interior, el resto conserva la geometría inicial

### target
```python
lumi.target(src, dst, target, target_func='abs2', error_func='abs', weight=1, f=None, wl=None, Q0=None)
```
- src: puerto y modo de origen, p. ej. "o1@0-" (modo 0 entrante)
- dst: puerto y modo de destino, p. ej. "o2@0+" (modo 0 saliente), o modo de volumen de un monitor de caja, p. ej. "o2@0" (sin dirección): en ese caso abs2 es la energía almacenada en el modo por unidad de densidad espectral incidente en `src`, independiente de la forma del pulso
- target: valor objetivo
- target_func: se aplica al parámetro S antes de compararlo con el objetivo, p. ej. "abs2" para la potencia
- weight: peso de este objetivo en la función de pérdida global
- f o wl: frecuencia o longitud de onda; especifique una (no ambas)
- error_func: función de error `abs, identity, -, relu, relu-'
- Q0: con un `dst` de modo de volumen, su Q debido a pérdidas distintas del acoplamiento con `src`. El valor pasa a ser la fracción de potencia en estado estacionario desde `src` hacia el modo, 4 g g0 / (g + g0)^2, donde g = abs2 anterior es la tasa de acoplamiento y g0 = 2 pi f / Q0. El objetivo 1 corresponde al acoplamiento crítico

### optimizer
```python
lumi.optimizer(contrasts, rates, stoploss=0, startloss=1, lowloss=None, momentum=0.7, init_dA=0.03, ckptat=None, gradckptat=None, lbfgs=False)
```
Argumentos:

# - iters: número máximo de iteraciones de optimización
- stoploss: detiene la optimización si la pérdida cae por debajo de este valor en el último contraste
- contrasts: contraste binario (de 0 a 1) de cada etapa. <1 permite cambios de topología (nuevos agujeros o islas). =1 congela la topología, p. ej. optimización de forma
- rates: por etapa. Se pasa a la siguiente etapa cuando la mejor pérdida de las últimas 10 iteraciones mejora la mejor anterior a ellas en menos de este valor. Cada etapa se reanuda a partir del mejor diseño de la etapa anterior
- startloss: estimación inicial de la pérdida
- lowloss: pérdida a la que el contraste alcanza 1. Por defecto, 0.2 * stoploss + 0.8 * startloss
# - highiter: iteración en la que el contraste alcanza 1. None lo incrementa solo en función de la pérdida
- momentum: momento del optimizador
# - contrast: contraste binario inicial (de 0 a 1). Ligeramente <1 permite cambios de topología (nuevos agujeros o islas). =1 congela la topología, p. ej. optimización de forma. El contraste aumenta automáticamente de forma lineal hasta 1 a medida que la pérdida disminuye hasta `lowloss`
- init_dA: fracción del área de diseño modificada por el primer paso. Fija el tamaño de paso inicial
- ckptat: guarda un punto de control de la optimización cada _ iteraciones
- gradckptat: checkpointing de gradientes
- lbfgs: usa L-BFGS con búsqueda lineal con retroceso (máximo 3 intentos; después finaliza la etapa) en las etapas con contraste 1 en lugar de momento

### apply_design
```python
lumi.apply_design(c0, sol)
```
Reemplaza las regiones de diseño del componente `c0` por los diseños optimizados de la solución `sol`.

### design_from_gds
```python
lumi.design_from_gds(path, i=1)
```
Importa el `i`-ésimo diseño optimizado desde `path` como componente de gdsfactory.

## Utilidades
### straight
```python
lumi.straight(length=10, width: float = 1, layer=(1, 0))
```
Componente rectangular de gdsfactory con los puertos "o1" (oeste) y "o2" (este).

### remote_solve
```python
lumi.remote_solve(path, api=None, poll=3, fields=False, gzip=False, wait_down=1800)
```
Ejecuta la simulación de la carpeta `path` en un servidor de Luminescent. Sube las entradas, transmite el registro en tiempo real y descarga los resultados en `path`.
El servidor genera la animación; también se descargan las instantáneas de campo (cortes de visualización).
Los planos de campo de los monitores de puerto (grandes: cada componente a cada longitud de onda) solo se descargan si `fields`.
Las subidas y descargas solo se comprimen con gzip si `gzip`.
Usa el token guardado por `set_token`.
`api` toma por defecto el servidor de Lumi más cercano.
Reutiliza una solución anterior con el mismo `watermark` (en esta carpeta, o de un trabajo suyo anterior en el servidor) en lugar de volver a ejecutarla.
También se guarda una copia de la ejecución en su directorio personal de Lumi Colab en ~/runs/<folder name>, que puede verse en la GUI de ejecuciones (Runs en la página de inicio de Lumi Colab).
Si el servidor se cae mientras se ejecuta el trabajo, sigue esperando hasta `wait_down` segundos a que vuelva (el trabajo se reanuda).
Lanza un error si hay una versión más reciente de luminescent en PyPI.

### load_sparams
```python
lumi.load_sparams(sparams)
```
Convierte diccionarios {"re", "im"} en valores complejos.

### version
```python
lumi.version()
```
Versión instalada del paquete de Python luminescent. Su X.Y debe coincidir con el X.Y.0 del binario del solver.

