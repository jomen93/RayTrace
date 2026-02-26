# C++ Ray Tracer Module

Módulo C++ acelerado para ray tracing de agujeros negros. ~10-20x más rápido que Python puro.

## Requisitos

- Python 3.6+
- pybind11: `pip install pybind11`
- Compilador C++14 (g++ o clang++)
- OpenMP (opcional pero recomendado): `brew install libomp` (macOS) o `apt-get install libgomp1` (Linux)

## Compilación

```bash
cd cpp
make

# Copiar módulo al directorio principal
cp raytracer_cpp.so ..
cd ..
```

## Uso

```bash
python main_cpp.py
```

## Solución de problemas

### macOS: "libomp not found"
```bash
brew install libomp
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
```

### Linux: "pybind11 not found"
```bash
pip install pybind11
```

### Error de compilación
```bash
make clean
make check  # Ver configuración
make
```

## Benchmarks esperados

| Resolución | Python | C++ | Speedup |
|------------|--------|-----|---------|
| 128×128 | ~3 min | ~10 sec | 18x |
| 256×256 | ~10 min | ~30 sec | 20x |
| 512×512 | ~41 min | ~2-4 min | 15x |
| 1024×1024 | ~3 horas | ~10 min | 18x |

## Arquitectura

- `raytracer.h` - Header con definiciones
- `raytracer.cpp` - Implementación del integrador RK4
- `binding.cpp` - Bindings Python con pybind11
- Integración numérica en C++ puro (~50x más rápido que scipy.integrate)
- Paralelismo con OpenMP
