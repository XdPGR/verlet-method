import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import sys
import matplotlib

matplotlib.use('TkAgg')

def create_animation(data_folder):
    data_folder = data_folder.rstrip('/')
    
    files = sorted(glob.glob(os.path.join(data_folder, "conf_*.csv")),
                   key=lambda x: int(x.split('_')[-1].split('.')[0]))
    
    if not files:
        print(f"Ошибка: нет файлов в {data_folder}")
        return
    
    print(f"Найдено {len(files)} файлов")
    
    first_data = np.loadtxt(files[0], delimiter=',', skiprows=1)
    if len(first_data.shape) == 1:
        first_data = first_data.reshape(1, -1)
    
    L = max(first_data[:, 0].max(), first_data[:, 1].max()) + 5
    
    fig, ax = plt.subplots(figsize=(7, 7))
    
    # Устанавливаем позицию окна
    try:
        mngr = plt.get_current_fig_manager()
        mngr.window.wm_geometry("+50+50")
    except:
        try:
            mngr.set_window_title("Анимация частиц")
            mngr.resize(700, 700)
            mngr.move(50, 50)
        except:
            pass
    
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_xlabel("X", fontsize=12)
    ax.set_ylabel("Y", fontsize=12)
    ax.set_title(f"Молекулярная динамика: {os.path.basename(data_folder)}", fontsize=14)
    ax.grid(True, alpha=0.3)
    
    scat = ax.scatter([], [], s=50, c='blue', alpha=0.7)
    
    def init():
        scat.set_offsets(np.empty((0, 2)))
        return scat,
    
    def update(frame):
        data = np.loadtxt(files[frame], delimiter=',', skiprows=1)
        if len(data.shape) == 1:
            data = data.reshape(1, -1)
        
        if len(data) > 0:
            points = data[:, :2]
            scat.set_offsets(points)
            ax.set_title(f"Шаг: {frame * 100}, частиц: {len(points)}")
        return scat, 
    
    ani = animation.FuncAnimation(fig, update, frames=len(files),
                                  init_func=init, blit=True, interval=50)
    
    plt.show()
    return ani

if __name__ == "__main__":
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    else:
        folder = input("Введите путь к папке с данными: ")
    create_animation(folder)
