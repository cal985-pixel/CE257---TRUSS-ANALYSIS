import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math

class Joint:
    def __init__(self, ID, x, y):
        self.ID = ID
        self.x = float(x)
        self.y = float(y)
        self.Px = 0.0      # Horizontal force (kN)
        self.Py = 0.0      # Vertical forcce (kN)
        self.support_type = None  # 'pin' or 'roller'

class Member:
    def __init__(self, ID, start, end, E=200e9, A=0.01):
        self.ID = ID
        self.joint_start = start
        self.joint_end = end
        self.E = float(E)   # Young’s Modulus(Pa)
        self.A = float(A)   # Area(m²)
        self.force = 0.0    # kN (positive = tension)

    def length(self):
        return math.hypot(self.joint_end.x - self.joint_start.x,
                          self.joint_end.y - self.joint_start.y)  

    def cos_angle(self):
        return (self.joint_end.x - self.joint_start.x) / self.length() #horizontal component / hypoteneuse

    def sin_angle(self):
        return (self.joint_end.y - self.joint_start.y) / self.length()  #verical component / hypoteneuse

#METHOD OF JOINTS SOLVER
class MethodOfJointsSolver:
    def __init__(self, joints, members):
        self.joints = joints
        self.members = members
        self.member_force = {}   # key: member object, value: force (kN)
        self.reaction_x = {}     # key: joint ID, value: reaction (kN)
        self.reaction_y = {}     # key: joint ID, value: reaction (kN)

    def solve(self):
        #Initialise unknowns to None 
        for m in self.members:
            self.member_force[m] = None
        for j in self.joints:
            self.reaction_x[j.ID] = None
            self.reaction_y[j.ID] = None

        mem_f = self.member_force.copy()  # Making copies to make repeated solving easier
        rx = self.reaction_x.copy()
        ry = self.reaction_y.copy()

        max_iterations = 100   #Solving for unknowns by repeatedly going through all joints
        for iteration in range(max_iterations):
            solved_anything = False

            for j in self.joints:  # Loop over every joint
                sum_Fx = 0.0  # Equilibrium equations: sum Fx = 0 , sum Fy = 0
                sum_Fy = 0.0

                sum_Fx = sum_Fx - j.Px  # Start with applied loads (negative because they go to RHS)
                sum_Fy = sum_Fy - j.Py

                # Collect unknowns at this joint
                unknown_list = []  # each item: (type, ID, coeff_x, coeff_y)
                # type can be 'M' for member, 'RX' for reaction x, 'RY' for reaction y

                for m in self.members:   # Add member forces that act on this joint
                    if m.joint_start is j:   # Member pulls from start to end
                        fx = m.cos_angle()
                        fy = m.sin_angle()
                        if mem_f[m] is None:
                            unknown_list.append(('M', m, fx, fy))
                        else:
                            sum_Fx = sum_Fx + mem_f[m] * fx
                            sum_Fy = sum_Fy + mem_f[m] * fy
                    # If this joint is the end of the member
                    elif m.joint_end is j:
                        fx = -m.cos_angle()    # Member pulls from end to start (opposite direction)
                        fy = -m.sin_angle()
                        if mem_f[m] is None:
                            unknown_list.append(('M', m, fx, fy))
                        else:
                            sum_Fx = sum_Fx + mem_f[m] * fx
                            sum_Fy = sum_Fy + mem_f[m] * fy

                if j.support_type == 'pin':      # Add support reactions at this joint
                    if rx[j.ID] is None:    # Horizontal reaction
                        unknown_list.append(('RX', j.ID, 1.0, 0.0))
                    else:
                        sum_Fx = sum_Fx + rx[j.ID]
                    if ry[j.ID] is None:     # Vertical reaction
                        unknown_list.append(('RY', j.ID, 0.0, 1.0))
                    else:
                        sum_Fy = sum_Fy + ry[j.ID]
                elif j.support_type == 'roller':
                    if ry[j.ID] is None:   # Only vertical reaction
                        unknown_list.append(('RY', j.ID, 0.0, 1.0))
                    else:
                        sum_Fy = sum_Fy + ry[j.ID]

                #we have a list of unknowns at this joint. solve if there are 1 or 2 unknowns
                if len(unknown_list) == 1:
                    # One unknown: use one equation (the one with non-zero coefficient)
                    typ, id1, cx, cy = unknown_list[0]
                    if abs(cx) > 0.000000001:   # Use x-equation
                        force = -sum_Fx / cx
                    elif abs(cy) > 0.000000001: # Use y-equation
                        force = -sum_Fy / cy
                    else:
                        raise ValueError("Joint has unknown with no coefficient!")   # Store the solved force
                    if typ == 'M':
                        mem_f[id1] = force
                    elif typ == 'RX':
                        rx[id1] = force
                    elif typ == 'RY':
                        ry[id1] = force
                    solved_anything = True

                elif len(unknown_list) == 2:   # Two unknowns: solve 2 equations
                    (typ1, id1, a1x, a1y) = unknown_list[0]
                    (typ2, id2, a2x, a2y) = unknown_list[1]

                    det = a1x * a2y - a2x * a1y    # Solve using Cramer's rule
                    if abs(det) < 0.000000001:   # Cannot solve yet (parallel members)
                        continue

                    F1 = (-sum_Fx * a2y + sum_Fy * a2x) / det
                    F2 = (-a1x * sum_Fy + a1y * sum_Fx) / det

                    # Store results
                    if typ1 == 'M':
                        mem_f[id1] = F1
                    elif typ1 == 'RX':
                        rx[id1] = F1
                    elif typ1 == 'RY':
                        ry[id1] = F1

                    if typ2 == 'M':
                        mem_f[id2] = F2
                    elif typ2 == 'RX':
                        rx[id2] = F2
                    elif typ2 == 'RY':
                        ry[id2] = F2

                    solved_anything = True

            if not solved_anything:   # If solve any new unknown in this full pass, stop
                break

        for m in self.members:  # After all iterations, check that every unknown is solved
            if mem_f[m] is None:
                raise ValueError("Not all member forces could be solved – truss may be indeterminate or unstable")
        for j in self.joints:
            if j.support_type == 'pin':
                if rx[j.ID] is None or ry[j.ID] is None:
                    raise ValueError("Reaction at pin support not solved")
            elif j.support_type == 'roller':
                if ry[j.ID] is None:
                    raise ValueError("Reaction at roller support not solved")

        self.member_force = mem_f    # Copy solved forces back to original storage
        self.reaction_x = rx
        self.reaction_y = ry

        for m in self.members:  #update the member objects themselves
            m.force = self.member_force[m]

        reactions = {}   # Build a reactions dictionary
        for j in self.joints:
            if j.support_type == 'pin':
                reactions[(j.ID, 'x')] = self.reaction_x[j.ID]
                reactions[(j.ID, 'y')] = self.reaction_y[j.ID]
            elif j.support_type == 'roller':
                reactions[(j.ID, 'y')] = self.reaction_y[j.ID]
        return reactions

    def compute_deflections(self): # Deflection calculation
        real_elong = {}   #ompute the real elongation of each member (in meters)
        for m in self.members:
            force_newtons = m.force * 1000.0   # convert kN to N
            if abs(force_newtons) < 0.000001:
                real_elong[m] = 0.0
            else:
                L = m.length()
                E = m.E
                A = m.A
                elongation = (force_newtons * L) / (E * A)   # meters
                real_elong[m] = elongation

        for j in self.joints:   #storing deflections in mm directly on each joint (as j.dx, j.dy) For each joint, compute horizontal deflection
            old_Px = {}   # Save the current loads (so we can restore later)
            old_Py = {}
            for jj in self.joints:
                old_Px[jj] = jj.Px
                old_Py[jj] = jj.Py

            for jj in self.joints:  # Apply virtual load: 1 kN in +x direction at joint j, zero elsewhere
                jj.Px = 0.0
                jj.Py = 0.0
            j.Px = 1.0

            virt_solver = MethodOfJointsSolver(self.joints, self.members)  # Solve the virtual truss
            virt_solver.solve()

            total = 0.0      # Sum (virtual force * real elongation) for all members
            for m in self.members:
                virt_force = virt_solver.member_force[m]   # in kN
                total = total + virt_force * real_elong[m]

            # Convert total (kN·m) to deflection in mm
            # 1 kN·m = 1000 N·m, and we want mm = m * 1000
            # Final formula: deflection_mm = total
            # (See derivation: total in kN·m -> divide by 1000 to get meters -> multiply by 1000 to get mm = total)
            j.dx = total   # in mm

            for jj in self.joints:  # Restore original loads
                jj.Px = old_Px[jj]
                jj.Py = old_Py[jj]

        for j in self.joints:  #compute vertical deflections (1 kN downward = +y)
            # Save loads
            old_Px = {}
            old_Py = {}
            for jj in self.joints:
                old_Px[jj] = jj.Px
                old_Py[jj] = jj.Py

            for jj in self.joints:   # Apply virtual load: 1 kN in +y at joint j
                jj.Px = 0.0
                jj.Py = 0.0
            j.Py = 1.0

            virt_solver = MethodOfJointsSolver(self.joints, self.members)
            virt_solver.solve()

            total = 0.0
            for m in self.members:
                virt_force = virt_solver.member_force[m]
                total = total + virt_force * real_elong[m]

            j.dy = total   # in mm

            for jj in self.joints:   # Restore loads
                jj.Px = old_Px[jj]
                jj.Py = old_Py[jj]

class InputDialog(tk.Toplevel):
    def __init__(self, parent, title, fields):
        super().__init__(parent)
        self.title(title)
        self.configure(bg="#16213e")
        self.resizable(False, False)
        self.result = None
        self.entries = {}
        tk.Label(self, text=title, bg="#16213e", fg="#e94560",
                 font=("Courier", 11, "bold")).pack(pady=(15, 10), padx=20)
        for field in fields:
            row = tk.Frame(self, bg="#16213e")
            row.pack(fill=tk.X, padx=20, pady=4)
            tk.Label(row, text=field + ":", bg="#16213e", fg="white",
                     font=("Courier", 9), width=14, anchor='w').pack(side=tk.LEFT)
            e = tk.Entry(row, bg="#0f3460", fg="white", insertbackground="white",
                         font=("Courier", 9), width=12)
            e.pack(side=tk.LEFT, padx=5)
            self.entries[field] = e
        tk.Button(self, text="OK", command=self.on_ok, bg="#e94560", fg="white",
                  font=("Courier", 9, "bold"), relief=tk.FLAT, width=10).pack(pady=15)
        self.grab_set()
        self.wait_window()
    def on_ok(self):
        self.result = {k: v.get() for k, v in self.entries.items()}
        self.destroy()

def ask_fields(parent, title, fields):
    d = InputDialog(parent, title, fields)
    return d.result

class TrussApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Truss Analysis - Method of Joints")
        self.root.configure(bg="#1a1a2e")
        self.root.geometry("1200x750")
        self.joints, self.members = [], []
        self.jc, self.mc = 1, 1
        self.grid_spacing = 1.0
        self.mode = "joint"
        self.mstart = None
        self.build_ui()

    def build_ui(self):
        left = tk.Frame(self.root, bg="#16213e", width=280)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(10, 0), pady=10)
        left.pack_propagate(False)
        tk.Label(left, text="TRUSS ANALYZER", bg="#16213e", fg="#e94560",
                 font=("Courier", 14, "bold")).pack(pady=(15, 5))
        tk.Label(left, text="Method of Joints", bg="#16213e", fg="#a8a8b3",
                 font=("Courier", 9)).pack(pady=(0, 15))
        gf = tk.LabelFrame(left, text="Grid", bg="#16213e", fg="#e94560", font=("Courier", 9, "bold"))
        gf.pack(fill=tk.X, padx=10, pady=5)
        tk.Label(gf, text="Spacing (m):", bg="#16213e", fg="white", font=("Courier", 9)).pack(side=tk.LEFT, padx=5)
        self.gv = tk.StringVar(value="1")
        tk.Entry(gf, textvariable=self.gv, width=6, bg="#0f3460", fg="white",
                 insertbackground="white").pack(side=tk.LEFT, padx=5)
        tk.Button(gf, text="Set", command=self.set_grid, bg="#e94560", fg="white",
                  font=("Courier", 8), relief=tk.FLAT).pack(side=tk.LEFT, padx=5, pady=5)
        mf = tk.LabelFrame(left, text="Mode", bg="#16213e", fg="#e94560", font=("Courier", 9, "bold"))
        mf.pack(fill=tk.X, padx=10, pady=5)
        for label, mode in [("Add Joint", "joint"), ("Add Member", "member"),
                             ("Apply Load", "load"), ("Add Support", "support")]:
            tk.Button(mf, text=label, command=lambda m=mode: self.set_mode(m),
                      bg="#0f3460", fg="white", font=("Courier", 9),
                      relief=tk.FLAT, width=18).pack(pady=3, padx=8)
        af = tk.LabelFrame(left, text="Actions", bg="#16213e", fg="#e94560", font=("Courier", 9, "bold"))
        af.pack(fill=tk.X, padx=10, pady=5)
        tk.Button(af, text="SOLVE", command=self.solve, bg="#e94560", fg="white",
                  font=("Courier", 10, "bold"), relief=tk.FLAT, width=18).pack(pady=5, padx=8)
        tk.Button(af, text="Show Deflection", command=self.show_deflection, bg="#0f3460",
                  fg="white", font=("Courier", 9), relief=tk.FLAT, width=18).pack(pady=3, padx=8)
        tk.Button(af, text="Clear All", command=self.clear_all, bg="#333",
                  fg="#a8a8b3", font=("Courier", 9), relief=tk.FLAT, width=18).pack(pady=3, padx=8)
        self.sv = tk.StringVar(value="Mode: Add Joint")
        tk.Label(left, textvariable=self.sv, bg="#16213e", fg="#e94560",
                 font=("Courier", 9), wraplength=250).pack(pady=10, padx=10)
        rf = tk.LabelFrame(left, text="Results", bg="#16213e", fg="#e94560", font=("Courier", 9, "bold"))
        rf.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        self.rt = tk.Text(rf, bg="#0a0a1a", fg="#00ff88", font=("Courier", 8),
                          height=10, state=tk.DISABLED)
        self.rt.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        right = tk.Frame(self.root, bg="#1a1a2e")
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.fig.patch.set_facecolor("#1a1a2e")
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.redraw()

    def set_grid(self):
        try:
            self.grid_spacing = float(self.gv.get())
            self.redraw()
        except ValueError:
            messagebox.showerror("Error", "Invalid grid spacing")

    def set_mode(self, mode):
        self.mode = mode
        self.mstart = None
        self.sv.set(f"Mode: {mode.capitalize()}")

    def on_click(self, event):
        if event.inaxes != self.ax:
            return
        x = round(event.xdata / self.grid_spacing) * self.grid_spacing
        y = round(event.ydata / self.grid_spacing) * self.grid_spacing
        actions = {"joint": self.add_joint, "member": self.add_member,
                   "load": self.add_load, "support": self.add_support}
        actions[self.mode](x, y)

    def snap(self, x, y):
        return next((j for j in self.joints if abs(j.x - x) <= 0.5 and abs(j.y - y) <= 0.5), None)

    def add_joint(self, x, y):
        if self.snap(x, y):
            self.sv.set("Joint already exists!")
            return
        data = ask_fields(self.root, "Add Joint", ["Joint Name"])
        if not data:
            return
        name = data["Joint Name"].strip() or f"J{self.jc}"
        j = Joint(name, x, y)
        self.jc += 1
        self.joints.append(j)
        self.redraw()
        self.sv.set(f"Added {j.ID} at ({x}, {y})")

    def add_member(self, x, y):
        j = self.snap(x, y)
        if not j:
            self.sv.set("Click on an existing joint!")
            return
        if not self.mstart:
            self.mstart = j
            self.sv.set(f"Start: {j.ID}. Now click end joint.")
        else:
            if self.mstart.ID == j.ID:
                self.sv.set("Cannot connect joint to itself!")
                self.mstart = None
                return
            for m in self.members:
                if (m.joint_start.ID == self.mstart.ID and m.joint_end.ID == j.ID) or \
                   (m.joint_start.ID == j.ID and m.joint_end.ID == self.mstart.ID):
                    self.sv.set("Member already exists!")
                    self.mstart = None
                    return
            data = ask_fields(self.root, "Add Member", ["Member Name", "E (Pa)", "A (m2)"])
            if not data:
                self.mstart = None
                return
            name = data["Member Name"].strip() or f"M{self.mc}"
            E = float(data["E (Pa)"] or 200e9)
            A = float(data["A (m2)"] or 0.01)
            mem = Member(name, self.mstart, j, E, A)
            self.mc += 1
            self.members.append(mem)
            self.sv.set(f"Added {mem.ID}: {self.mstart.ID} to {j.ID}")
            self.mstart = None
            self.redraw()

    def add_load(self, x, y):
        j = self.snap(x, y)
        if not j:
            self.sv.set("Click on an existing joint!")
            return
        data = ask_fields(self.root, f"Load at {j.ID}", ["Px (kN)", "Py (kN)"])
        if data:
            j.Px = float(data["Px (kN)"] or 0)
            j.Py = float(data["Py (kN)"] or 0)
            self.sv.set(f"Load at {j.ID}: ({j.Px}, {j.Py}) kN")
            self.redraw()

    def add_support(self, x, y):
        j = self.snap(x, y)
        if not j:
            self.sv.set("Click on an existing joint!")
            return
        data = ask_fields(self.root, f"Support at {j.ID}", ["Type (pin/roller)"])
        if data and data["Type (pin/roller)"] in ['pin', 'roller']:
            j.support_type = data["Type (pin/roller)"]
            self.sv.set(f"{j.support_type} at {j.ID}")
            self.redraw()
        else:
            self.sv.set("Use pin or roller only!")

    def redraw(self, force_colors=False):
        self.ax.clear()
        self.ax.set_facecolor("#0f3460")
        for i in np.arange(-10, 10 + self.grid_spacing, self.grid_spacing):
            self.ax.axhline(i, color="#1a3a6e", linewidth=0.5)
            self.ax.axvline(i, color="#1a3a6e", linewidth=0.5)
        self.ax.set_xlim(-1, 10)
        self.ax.set_ylim(-1, 8)
        self.ax.set_aspect('equal')
        self.ax.tick_params(colors='white')
        for sp in ['bottom', 'left']:
            self.ax.spines[sp].set_color('#e94560')
        for sp in ['top', 'right']:
            self.ax.spines[sp].set_color('#0f3460')
        for mem in self.members:
            color = ("#ff4444" if mem.force < 0 else "#44ff88") if force_colors else "#00aaff"
            self.ax.plot([mem.joint_start.x, mem.joint_end.x],
                         [mem.joint_start.y, mem.joint_end.y], color=color, linewidth=2.5, zorder=2)
            mx = (mem.joint_start.x + mem.joint_end.x) / 2
            my = (mem.joint_start.y + mem.joint_end.y) / 2
            lbl = f"{abs(mem.force):.1f}kN({'C' if mem.force<0 else 'T'})" if force_colors else mem.ID
            self.ax.text(mx, my + 0.15, lbl, color=color, fontsize=7, ha='center', fontweight='bold')
        for j in self.joints:
            self.ax.plot(j.x, j.y, 'o', color="#e94560" if j.support_type else "#00ff88",
                         markersize=8, zorder=4)
            self.ax.text(j.x + 0.1, j.y + 0.2, j.ID, color="white", fontsize=8, fontweight='bold')
            if j.support_type == 'pin':
                self.ax.plot(j.x, j.y - 0.3, '^', color="#ffaa00", markersize=12)
                self.ax.text(j.x, j.y - 0.55, 'PIN', color="#ffaa00", fontsize=6, ha='center')
            elif j.support_type == 'roller':
                self.ax.plot(j.x, j.y - 0.3, 'o', color="#ffaa00", markersize=10,
                             markerfacecolor='none', markeredgewidth=2)
                self.ax.text(j.x, j.y - 0.55, 'ROLLER', color="#ffaa00", fontsize=6, ha='center')
            if j.Px != 0 or j.Py != 0:
                # arrow direction
                mag = math.hypot(j.Px, j.Py) + 0.001
                dx_arrow = j.Px / mag * 0.5
                dy_arrow = j.Py / mag * 0.5
                self.ax.annotate('', xy=(j.x + dx_arrow, j.y + dy_arrow),
                                 xytext=(j.x, j.y),
                                 arrowprops=dict(arrowstyle='->', color='#ff6600', lw=2))
                self.ax.text(j.x+0.15, j.y-0.4,
                             f"({j.Px:.1f},{j.Py:.1f})kN", color="#ff6600", fontsize=6)
        self.canvas.draw()

    def solve(self):
        if not self.joints or not self.members:
            messagebox.showerror("Error", "Add joints and members first!")
            return
        try:
            solver = MethodOfJointsSolver(self.joints, self.members)
            reactions = solver.solve()
            out = "MEMBER FORCES (kN)\n"
            for mem in self.members:
                out += f"{mem.ID}: {mem.force:+.3f} kN ({'TENSION' if mem.force > 0 else 'COMPRESSION'})\n"
            out += "\nREACTIONS (kN)\n"
            for (jid, dir), val in reactions.items():
                out += f"{jid} ({dir}): {val:+.3f} kN\n"
            self.rt.config(state=tk.NORMAL)
            self.rt.delete(1.0, tk.END)
            self.rt.insert(tk.END, out)
            self.rt.config(state=tk.DISABLED)
            self.sv.set("Solved successfully (Method of Joints)!")
            self.redraw(force_colors=True)
        except Exception as e:
            messagebox.showerror("Solver Error", str(e))

    def show_deflection(self):
        if not any(m.force != 0 for m in self.members):
            messagebox.showwarning("Warning", "Solve the truss first!")
            return
        # Compute deflections using unit load method
        try:
            solver = MethodOfJointsSolver(self.joints, self.members)
            solver.solve()  # ensure forces are up to date
            solver.compute_deflections()
            fig2, ax2 = plt.subplots(figsize=(8, 6))  # Now plot deformed shape
            fig2.patch.set_facecolor("#1a1a2e")
            ax2.set_facecolor("#0f3460")
            ax2.set_title("Deformed Shape (Method of Joints + Unit Load)", color="white", fontsize=12)
            scale = 200  # magnification factor for deflections (mm to m)
            # Draw original members
            for mem in self.members:
                ax2.plot([mem.joint_start.x, mem.joint_end.x],
                         [mem.joint_start.y, mem.joint_end.y], '--', color="#444", linewidth=1)
            # Draw deformed members
            for mem in self.members:
                xs = [mem.joint_start.x + mem.joint_start.dx / 1000 * scale,
                      mem.joint_end.x + mem.joint_end.dx / 1000 * scale]
                ys = [mem.joint_start.y + mem.joint_start.dy / 1000 * scale,
                      mem.joint_end.y + mem.joint_end.dy / 1000 * scale]
                ax2.plot(xs, ys, '-', color="#e94560", linewidth=2.5)
            # Plot joints
            for j in self.joints:
                ax2.plot(j.x, j.y, 'bo', markersize=6, label='Original' if j==self.joints[0] else "")
                ax2.plot(j.x + j.dx/1000*scale, j.y + j.dy/1000*scale,
                         'ro', markersize=5, label='Deformed' if j==self.joints[0] else "")
                ax2.text(j.x, j.y, j.ID, color='white', fontsize=8)
            ax2.tick_params(colors='white')
            ax2.set_aspect('equal')
            ax2.legend(["Original", "Deformed"], facecolor="#16213e", labelcolor="white")
            plt.tight_layout()
            plt.show()
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def clear_all(self):
        self.joints.clear()
        self.members.clear()
        self.jc = self.mc = 1
        self.mstart = None
        self.rt.config(state=tk.NORMAL)
        self.rt.delete(1.0, tk.END)
        self.rt.config(state=tk.DISABLED)
        self.sv.set("Cleared. Mode: Add Joint")
        self.redraw()

if __name__ == "__main__":
    root = tk.Tk()
    app = TrussApp(root)
    root.mainloop()