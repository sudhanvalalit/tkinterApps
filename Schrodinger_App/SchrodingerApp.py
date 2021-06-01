import tkinter as tk
from tkinter import Scrollbar, Frame
from tkinter.tix import *
from tkinter.constants import BOTH, RIGHT, VERTICAL
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import seaborn as sns
import os
import shutil

sns.set()

hbar = 1.0  # 197.33e3
m = 1.0  # 511.0
N = 100
Nmax = 1
epsilon = 1e-10
energyScale = 1.0  # 3.5
xmin, xmax = 0.0, 6.0

# Main window
window = tk.Tk()
window.title("Schrodinger Equation Solver")
window.geometry('1000x800+700+200')
# window.pack()

# gettig potential and l value as input
pot_var = tk.StringVar()
Nmax_var = tk.StringVar()
xmin_var = tk.StringVar()
xmax_var = tk.StringVar()


def schro_var():
    '''
    Gets variables from the textbox on the screen.
    '''
    potential = pot_var.get()
    global Nmax
    Nmax = int(Nmax_var.get())
    global xmin, xmax
    xmin, xmax = float(xmin_var.get()), float(xmax_var.get())
    return potential


def V(x):
    '''
    Potential function. Reads potential from the screen and evaluates on the fly.
    '''
    potential = schro_var()
    potential = eval(potential)
    return potential


def buttonPushed():
    '''
    The exit button function.
    '''
    window.destroy()


def plot2d(x, data, my_canvas, eigenEnergies):
    '''
    Makes plots of the wavefunctions and displays corresponding eigenvalues on screen.
    '''
    N = len(data)
    figSize = 3.5
    rowNumber = 0
    dest = 'Schrodinger_App/plots/'
    if not os.path.exists(dest):
        os.makedirs(dest)
    for i in range(N):
        if i % 4 == 0:
            rowNumber += 1
        fig = plt.figure(figsize=(figSize, figSize), dpi=100)
        fileName = f"Wavefunction_{i+1}.png"
        plotTitle = f"Eigenvalue {i+1}: {eigenEnergies[i]:.2f}"
        plt.plot(x, data[i])
        plt.savefig(fileName)
        plt.title(plotTitle)
        canvas = FigureCanvasTkAgg(fig, master=my_canvas)
        canvas.draw()
        canvas.get_tk_widget().grid(row=rowNumber, column=i % 4)
        shutil.move(fileName, dest)

    return canvas


def eval_eigenvalues():
    '''
    This method of evaluating eigenvalues is called Numerov method. The method used 
    below can be found in N. Zettili -- Quantum Mechanics.
    '''
    potential = schro_var()
    psi = np.zeros(N)
    psi_final = []
    eigenEnergies = np.zeros(Nmax)
    psi_left, psi_right = 1.0e-6, 0.0
    psi[0], psi[1] = psi_left, psi_left + 1.0e-6
    ElowerLimit = 0.0
    EupperLimit = 10.0
    h0 = (xmax - xmin) / N
    xRange = np.arange(xmin, xmax, h0)
    Endsign = -1
    for nQuantum in range(Nmax):
        Limits_are_defined = False
        while Limits_are_defined is False:
            nodes_plus = 0
            Etrial = EupperLimit
            for i in range(2, N):
                Ksquare = 2.0 * (Etrial - V(xRange[i]))
                psi[i] = (
                    2.0
                    * psi[i - 1]
                    * (1.0 - 5.0 * h0 * h0 * Ksquare / 12.0)
                    / (1.0 + (h0 * h0 * Ksquare / 12.0))
                    - psi[i - 2]
                )
                if psi[i] * psi[i - 1] < 0:
                    nodes_plus += 1
            if EupperLimit < ElowerLimit:
                EupperLimit = np.max([2 * EupperLimit, -2.0 * EupperLimit])
            if nodes_plus > nQuantum + 1:
                EupperLimit *= 0.7
            elif nodes_plus < nQuantum + 1:
                EupperLimit *= 2.0
            else:
                Limits_are_defined = True
        Endsign *= -1
        while EupperLimit - ElowerLimit > epsilon:
            Etrial = (EupperLimit + ElowerLimit) / 2.0
            for i in range(2, N):
                Ksquare = 2.0 * (Etrial - V(xRange[i]))
                psi[i] = (
                    2.0
                    * psi[i - 1]
                    * (1.0 - 5.0 * h0 * h0 * Ksquare / 12.0)
                    / (1.0 + (h0 * h0 * Ksquare / 12))
                    - psi[i - 2]
                )
            if Endsign * psi[-1] > psi_right:
                ElowerLimit = Etrial
            elif Endsign * psi[-1] < psi_right:
                EupperLimit = Etrial
            else:
                exit()

        Etrial = (EupperLimit + ElowerLimit) / 2
        eigenEnergies[nQuantum] = Etrial
        EupperLimit = Etrial
        ElowerLimit = Etrial

        # Normalization
        Integral = 0.0
        for i in range(N):
            Integral += 0.5 * h0 * (psi[i - 1] * psi[i - 1] + psi[i] * psi[i])

        normCoeff = np.sqrt(1.0 / Integral)
        psi_final.append(normCoeff * psi)

    potential = np.zeros(len(xRange))
    potential = [V(xRange[i]) for i in range(len(xRange))]
    # psi_final.append(potential)
    plot2d(xRange, psi_final, frame, eigenEnergies)


def reset_values():
    pot_var.set("")
    Nmax_var.set("")
    xmin_var.set("")
    xmax_var.set("")
    for widget in window.winfo_children():
        widget.destroy()

    main()


def main():
    global my_canvas, frame
    my_canvas = tk.Canvas(
        window, width=600, height=400)
    # my_canvas = FigureCanvasTkAgg(master=frame)
    my_canvas.pack(side=tk.TOP, expand=True, fill=BOTH)
    vbar = Scrollbar(my_canvas, orient=VERTICAL, command=my_canvas.yview)
    vbar.pack(side=RIGHT, fill=Y)
    # # frame.config(width=2000, height=2000)
    my_canvas.configure(yscrollcommand=vbar.set)
    my_canvas.bind("<Configure>", lambda e: my_canvas.configure(
        scrollregion=my_canvas.bbox("all")))
    frame = Frame(my_canvas, width=1000, height=800, bd=5)
    frame.pack(expand=True, fill=BOTH)
    my_canvas.create_window((0, 0), anchor='nw', window=frame)
    #
    potential_label = tk.Label(
        window, text="Potential function (in terms of x):", font=("helvetica", 13, "bold")
    )
    Nmax_label = tk.Label(
        window, text="No. of Eigenvalues:", font=("helvetica", 13, "bold"))
    xmin_label = tk.Label(window, text="xmin:", font=("helvetica", 13, "bold"))
    xmax_label = tk.Label(window, text="xmax:", font=("helvetica", 13, "bold"))

    # creating entry for inpput

    potential_entry = tk.Entry(window, textvariable=pot_var)
    Nmax_entry = tk.Entry(window, textvariable=Nmax_var)
    xmin_entry = tk.Entry(window, textvariable=xmin_var)
    xmax_entry = tk.Entry(window, textvariable=xmax_var)

    # submit button
    submit = tk.Button(text="Submit", command=lambda: eval_eigenvalues())

    # reset button

    reset = tk.Button(text="Reset", command=lambda: reset_values())

    # Create exit button

    exitButton = tk.Button(window, text="Exit",
                           relief=tk.RIDGE, command=buttonPushed)
    exitButton.pack(side=tk.RIGHT, padx=10)

    # placing of widgets
    potential_label.pack(side=tk.LEFT)
    potential_entry.pack(side=tk.LEFT, padx=10)
    Nmax_label.pack(side=tk.LEFT, padx=10)
    Nmax_entry.pack(side=tk.LEFT, padx=10)
    xmin_label.pack(side=tk.LEFT, padx=10)
    xmin_entry.pack(side=tk.LEFT, padx=10)
    xmax_label.pack(side=tk.LEFT, padx=10)
    xmax_entry.pack(side=tk.LEFT, padx=10)
    submit.pack(side=tk.LEFT, padx=10)
    reset.pack(side=tk.LEFT, padx=10)

    window.mainloop()


if __name__ == "__main__":
    main()
