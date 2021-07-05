import math

st = open("stars.txt", "r")
angles = []
cnt = 0
radius = math.atan(100 / 390000)
stars = 0
amount = 3597437 // 1

print(radius)

det = 500000
lst = [0 for i in range(det)]

while cnt < amount:
    cnt += 1
    ss = st.readline().strip()[:-2]
    star = list(map(float, ss.strip().split()))
    if len(star) != 4:
        break
    tg = (star[0] ** 2 + star[1] ** 2) ** 0.5 / star[2]
    ang = math.atan(tg)
    angles.append(ang)

    lst[int((2*ang/math.pi) * det)] += 1

    if ang < radius:
        stars += 1




from tkinter import *

root = Tk()
root.geometry('1500x900')

canv = Canvas(root, bg='white')
canv.pack(expand=True, fill=BOTH)


args = {'zero_x': 100, 'zero_y': 700, 'max_x': 1400, 'max_y': 100, 'mark_size': 5, 'r': 4, 'r_big': 10}
args['edge_x'] = args['max_x'] + 50
args['edge_y'] = args['max_y'] - 50

canv.create_line(args['zero_x'], args['zero_y'], args['edge_x'], args['zero_y'], width=3, arrow=LAST)  # Ox
canv.create_line(args['zero_x'], args['zero_y'], args['zero_x'], args['edge_y'], width=3, arrow=LAST)  # Oy
#canv.create_line(args['zero_x'], args['max_y'], args['edge_x'], args['max_y'], width=1, dash=(5, 4))  # y = 1
#canv.create_line(args['zero_x'], 2 * args['zero_y'] - args['max_y'], args['edge_x'],
  #                   2 * args['zero_y'] - args['max_y'], width=1, dash=(5, 4))  # y = -1
#canv.create_line(args['zero_x'], args['zero_y'], args['zero_x'],
 #                    2 * args['zero_y'] - args['edge_y'], width=3)  # also Oy

r = 1.3
mx = math.log(1 + max(lst))


for i in range(det // 2):
    x, y = args['zero_x'] + (args['edge_x'] - args['zero_x']) * math.log(1 + i) / math.log(1 + det), args['zero_y'] + (args['edge_y'] - args['zero_y']) * math.log(1 + lst[i]) / mx
    canv.create_oval(x - r, y - r, x + r, y + r, fill= 'black', outline = 'black')

for i in range(det // 2, det, 8):
    x, y = args['zero_x'] + (args['edge_x'] - args['zero_x']) * math.log(1 + i) / math.log(1 + det), args['zero_y'] + (args['edge_y'] - args['zero_y']) * math.log(1 + lst[i]) / mx
    canv.create_oval(x - r, y - r, x + r, y + r, fill= 'black', outline = 'black')


x_r = args['zero_x'] + (args['edge_x'] - args['zero_x']) * math.log(1 + int(math.pi * det / 600)) / math.log(1 + det)
canv.create_line(x_r, 0, x_r, 10000)

root.mainloop()

