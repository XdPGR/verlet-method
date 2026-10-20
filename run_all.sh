#!/bin/bash

cd "$(dirname "$0")"

echo "=========================================="
echo "    МОЛЕКУЛЯРНАЯ ДИНАМИКА"
echo "=========================================="
echo ""
echo "Выберите состояние:"
echo "  1 - ГАЗ (случайные положения)"
echo "  2 - ГАЗ (регулярная решетка)"
echo "  3 - ЖИДКОСТЬ (малая)"
echo "  4 - ЖИДКОСТЬ (большая)"
echo "  5 - ТВЕРДОЕ ТЕЛО (квадратная решетка)"
echo "  6 - ТВЕРДОЕ ТЕЛО (треугольная решетка)"
echo ""
read -p "Введите номер (1-6): " choice

SOURCE_FILE="src/main_verle.c"

if [ ! -f "$SOURCE_FILE" ]; then
    echo "ОШИБКА: Не найден $SOURCE_FILE"
    exit 1
fi

case $choice in
    1) STATE="STATE_GAS"; SIM_NAME="sim_gas"; DATA_FOLDER="data/gas"; STATE_NAME="ГАЗ (случайные положения)";;
    2) STATE="STATE_GAS_REGULAR"; SIM_NAME="sim_gas_regular"; DATA_FOLDER="data/gas_regular"; STATE_NAME="ГАЗ (регулярная решетка)";;
    3) STATE="STATE_LIQUID"; SIM_NAME="sim_liquid"; DATA_FOLDER="data/liquid"; STATE_NAME="ЖИДКОСТЬ (малая)";;
    4) STATE="STATE_LIQUID_LARGE"; SIM_NAME="sim_liquid_large"; DATA_FOLDER="data/liquid_large"; STATE_NAME="ЖИДКОСТЬ (большая)";;
    5) STATE="STATE_SOLID"; SIM_NAME="sim_solid"; DATA_FOLDER="data/solid"; STATE_NAME="ТВЕРДОЕ ТЕЛО (квадратная решетка)";;
    6) STATE="STATE_SOLID_TRI"; SIM_NAME="sim_solid_tri"; DATA_FOLDER="data/solid_tri"; STATE_NAME="ТВЕРДОЕ ТЕЛО (треугольная решетка)";;
    *) echo "Неверный выбор!"; exit 1;;
esac

mkdir -p "$DATA_FOLDER"

echo ""
echo "==> Компиляция..."
gcc -D$STATE -o "$SIM_NAME" "$SOURCE_FILE" -lm

if [ ! -f "$SIM_NAME" ]; then
    echo "Ошибка компиляции!"
    exit 1
fi

echo "==> Запуск симуляции: $STATE_NAME"
./"$SIM_NAME"

echo ""
echo "==> Запуск визуализации..."
python3 scripts/animate.py "$DATA_FOLDER" &
ANIM_PID=$!

echo ""
echo "Нажмите Enter, чтобы закрыть окно анимации и завершить программу..."
read

kill $ANIM_PID 2>/dev/null
pkill -f "animate.py" 2>/dev/null

echo ""
echo "Готово!"
