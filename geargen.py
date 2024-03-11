#
#  Generate DXF files gears
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  ---------------------------------------------------------------------------
#
#  Run 'python geargen.py --help' for instructions.
#

import argparse
import getopt, sys
from math import *

# internal gear
internal = 0

# number of teeth
N = 5

# module
m = 1.0

# clearance
c = 0.0

# addendum offset
off_a = 0.0

# addendum/dedendum multiplier
mul_a = 0.5
mul_d = 0.5

# pitch
p = 0.0

# spacing
s = 0.0

# tip relief
tr = 0.0

# x- and y-offset, rotation
off_x, off_y, off_phi = 0.0, 0.0, 0.0

# skip dxf header/footer generation
e_mode = 0

# layer name for generated dxf data
layername = ""

# output filename
outfile = ""

def check_positive(cast):
    def check(string):
        #just to please the IDE
        val=0
        try:
            val = cast(string)
        except ValueError:
            msg="%r needs to be a valid number!"
            raise argparse.ArgumentTypeError(msg)
        if val < 0:
            mst="%r is negative, but positive number is required!"
            raise argpars.ArgumentTypeError(msg)
        return val
    return check

parser = argparse.ArgumentParser(
    usage="usage: {0} [options]".format(sys.argv[0]),\
    epilog="Argument types as used in above options reference:\n\
    N = positive integer\n\
    U = unsigned decimal number\n\
    X = signed decimal number\n\
    T = text",
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", "--invert",
    dest="internal",
    action="store_true",
    help="generate outline for internal gear \
    (i.e. the meaning of pitch and clearance is inverted,\
    the sign of the spacing is changed and the tip relief\
    is reset to 0.0)")
parser.add_argument("-C", "--cycloid",
    dest="involute",
    action="store_false",
    help="generate cycloid gear instead of \
    involute. Does NOT accept all arguments!\n\
    wikipedia.org/wiki/Cycloid_gear")
parser.add_argument("-n", "--number",
    dest="N",
    default=12,
    metavar="N",
    type=check_positive(int),
    help="number of teeth (default = 5)")
parser.add_argument("-m", "--module",
    dest="m",
    default=2*pi,
    metavar="U",
    type=check_positive(float),
    help="module (default = 1.0)\n\
    (for gears to match, their modules must be equal. the gears\
    working diameter is [number of teeth] times [module])")
parser.add_argument("-c", "--clearance",
    dest="c",
    default=0.0,
    metavar="X",
    type=float,
    help="clearance (default = 0.0)\n\
    (move the gear bottom land down by this number of millimeters)")
parser.add_argument("-a", "--add-offset",
    dest="off_a",
    default=0.0,
    metavar="X",
    type=float,
    help="addendum offset (default = 0.0)\n\
    (move the gear top land up by this number of millimeters by\
    incrasing the addendum by this value. use with care!)")
parser.add_argument("-t", "--top",
    dest="t",
    default=0.0,
    metavar="X",
    type=check_positive(float),
    help="move the gear top land up by this number of millimeters by\
    linearly extending the top land. \
    WILL BREAK IF USED WITH INSUFFICIENT CLEARANCE! \
    (on corresponding wheel) \
    do not confuse with '-a'!")
parser.add_argument("-p", "--attack",
    dest="p",
    default=25.0,
    metavar="X",
    type=float,
    help="pitch (default = 20.0)\n\
    (Angle of attack at pitch-radius in degress.\
    Will be ignored if angle of attack in radians is set. Do not confuse with '-a'!)")
parser.add_argument("-P", "--Attack",
    dest="P",
    default=0.0,
    metavar="X",
    type=float,
    help="Pitch (no default)\n\
    (Angle of attack at pitch-radius in radians. Do not confuse with '-a'!)")
parser.add_argument("-s", "--spacing",
    dest="s",
    default=0.0,
    metavar="X",
    type=float,
    help="spacing (default = 0.0)\n\
    (move the flanks inside by this number of millimeters.)")
parser.add_argument("-T", "--relief",
    dest="tr",
    default=0.0,
    metavar="U",
    type=check_positive(float),
    help="tip relief (default = 0.0)\n\
    (cut off a triangle with this side length from the gear tooths)")
parser.add_argument("-D", "--dedendum",
    dest="mul_d",
    default=1.0,
    metavar="U",
    type=check_positive(float),
    help="dedendum multiplier (default = 0.5)")
parser.add_argument("-A", "--addendum",
    dest="mul_a",
    default=1.0,
    metavar="U",
    type=check_positive(float),
    help="addendum multiplier (default = 0.5)\n\
    (multiply module with this value to calculate dedendum/addendum)")
parser.add_argument("-e", "--headerless",
    dest="e_mode",
    action="store_true",
    help="generate only DXF entities without DXF header or DXF footer\
    (usefull for creating larger dxf files from a script)")
parser.add_argument("-l", "--layer",
    dest="layername",
    metavar="T",
    type=str,
    help="name of the DXF layer for generated entities\n\
    (without this option, no layer information is added to the DXF)")
parser.add_argument("-x",
    dest="off_x",
    default=0.0,
    metavar="X",
    type=float,
    help="add this value to all X-coordinates in the generated DXF file")
parser.add_argument("-y",
    dest="off_y",
    default=0.0,
    metavar="X",
    type=float,
    help="add this value to all Y-coordinates in the generated DXF file")
parser.add_argument("-r", "--rotate",
    dest="off_phi",
    default=0.0,
    metavar="X",
    type=float,
    help="rotate the generate wheel by this angle in degrees\n\
    (this can be used together with -e and -l to generate more complex\
    DXF files from scripts. default oriantation is gear at 0/0 with\
    a bottom land facing to the positive x-axis)")
parser.add_argument("-f", "--filename",
    dest="outfile",
    metavar="T",
    type=str,
    help="name of the DXF file to generate\n\
    (without this option, a (tk) window is opend and the generated gear\
    is displayed. note that the gear is scaled to fit the viewport.)")

args=parser.parse_args()

pitch=float(args.m*args.N)/(2*pi)

linedata = []

def line(x1, y1, x2, y2):
    global linedata
    linedata.append(x1)
    linedata.append(y1)
    linedata.append(x2)
    linedata.append(y2)

def dot(x,y):
    global linedata
    linedata.append(x)
    linedata.append(y)

def drop_line(n=1):
    global linedata
    linedata[-n*4:]=[]

def drop_dot(n=1):
    global linedata
    linedata[-n*2]=[]

def find_intersection(lines1, lines2):
    #SKIZZE!
    test1 = lines1[-1] - lines[0]
    test2 = lines2[-1] - lines[0]
    #if test1[0]==

#requires number of teeth, passed by args
def singleArc(phi_start, phi_stop, phi_step):
    ###constants
    #scale parameters to 0-1 range
    scale=phi_step/(phi_stop-phi_start)
    #pitch radius
    pitcher=args.N/(2*pi)

    ###variables
    psi=0
    invol_start=0
    steps=0
    before=()
    point=rack(0)
    after=rack(scale)
    x=y=0
    psi = []
    last=0
    if after[1]!=point[1]:
        invder=(after[0]-point[0])/(after[1]-point[1])
        #x=pitcher+point[0]; psi=(1-x/pitcher)*invder+...
        psi= [-point[0]/pitcher*invder + point[1]/pitcher]
        x=(point[0]+pitcher)*cos(psi[-1]) - point[0]*invder*sin(psi[-1])
        y=(point[0]+pitcher)*sin(psi[-1]) + point[0]*invder*cos(psi[-1])

    for i in range(1,int((phi_stop-phi_start)/phi_step)-1):
        before=point
        point=after
        after=rack((i+1)*scale)
        if before[1]==point[1] and after[1]==point[1]:
            steps +=1
            continue
        elif before[1]==point[1]:
            psi_stop,nil=involute_intersect_with_r2(invol_start, point[0])
            involute(invol_start, psi, psi_stop, steps)
            px,py=point
            steps=0
            invder=(after[0]-point[0])/(after[1]-point[1])
        elif after[1]==point[1]:
            invol_start=point[0]
            invder=(point[0]-before[0])/(point[1]-before[1])
        else:
            if before[0]==point[0] or point[0]==after[0]:
                invder=0
            else:
                inv=(point[1]-before[1])*(after[0]-point[0])+(after[1]-point[1])*(point[0]-before[0])
                if inv == 0:
                    invol_start=point[0]
                    continue
                invder=(2*(point[0]-before[0])*(after[0]-point[0])/inv)
        #x=pitcher+point[0]; psi=(1-x/pitcher)*invder+...
        nextpsi=(point[1]-point[0]*invder)/pitcher
        psi.append(nextpsi)
        print(point[1]-point[0]*invder)

        x=(pitcher+point[0])*cos(psi[-1]) - point[0]*invder*sin(psi[-1])
        y=(pitcher+point[0])*sin(psi[-1]) + point[0]*invder*cos(psi[-1])

        dot(x,y)
    print(pitcher)

#low = half angle of lower radius
#diag = half angle of tooth at max. width
#low*r == r/tan(diag)
#diag + low = pi/args.N - 1/r * sqrt((r/cos(diag))**2 - r**2)
#diag + low = pi/args.N - sqrt(1/cos(diag)**2 - 1)
#pi/args.N = diag+low - sqrt((1-cos(diag)**2)/cos(diag)**2)
#pi/args.N = diag+low - tan(diag) = diag + 1/tan(diag) - tan(diag)
#use low, requires just one atan. Using diag requires 2 tan-calls.
#pi/args.N = low + 1/atan(low) - 1/low

class Cog(object):

    def __init__(self, args):
        for k,v in args.__dict__.items():
            setattr(self, k, v)
        if args.internal:
            self.tip, self.clearance = args.c, args.t
            self.spacing = args.s * -1.0
        else:
            self.clearance, self.tip = args.c, args.t
            self.spacing = args.s
        self.pitch=args.m*args.N/(2*pi)

    def draw_wheel(self):
        cogSpan=2*pi/args.N
        for i in range(args.N):
            phi=cogSpan*i + args.off_phi*(pi/180)
            x,y=self.draw(phi, 100)

class Involute(Cog):

    def __init__(self, args):
        super(Involute, self).__init__(args)
        if self.P==0:
            self.P=args.p*pi/180
        self.D=args.mul_d
        self.A=args.mul_a
        self._variable_attack_involute_span()

    def _t_from_phi(self, inv, d=1e-8, runs=10000):
        top=pi/2
        low=-pi/2
        for i in range(runs):
            if top-low < d:
                return (top+low)/2
            elif (top+low)/2 - atan((top+low)/2) > inv:
                top=(top+low)/2
            else:
                low =(top+low)/2
        raise ValueError("Could not find suitably accurate alpha.")

    def _variable_attack_involute_span(self):
        """
        Calculates angle-span of involute, (=span)
        where to start drawing, (=shift)
        rotation angle of involute (=offset)
        """
        self.base=self.pitch*cos(self.P)
        self.bottom=self.pitch-self.D
        self.top=self.pitch+self.A
        if self.base > self.bottom:
            raise ValueError("Dedendum undercuts base circle by %r"%(self.base-self.bottom))
        maxHalfSpan=pi/self.N/2 #max angle btw. top a. bottom, one side of single tooth
        #all following are parameters of intersection with [varname]
        low = sqrt((self.bottom/self.base)**2 -1)
        top = sqrt((self.top   /self.base)**2 -1)
        pitch=sqrt((self.pitch /self.base)**2 -1)
        #spanned angle at given parameter
        ilow = low - atan(low)
        itop=  top - atan(top)
        ipitch=pitch-atan(pitch)

        if itop-ipitch > maxHalfSpan:
            upper=self._t_from_phi(ipitch+maxHalfSpan)
            self.top=sqrt(upper**2+1)*self.base
        else:
            upper=top
        if ipitch-ilow > maxHalfSpan:
            self.shift=self._t_from_phi(ipitch-maxHalfSpan)
            self.bottom=sqrt(self.shift**2+1)*self.base
        else:
            self.shift=low
        self.span=upper-self.shift
        self.offset=max(0,maxHalfSpan-(ipitch-ilow))
        self.low=ilow
        return self.span, self.shift, self.offset

    def draw(self, phi_start, steps):
        if self.tip!=0 and self.bottom == self.pitch-self.D:
            rt=self.top+self.tip
            frat=rt/self.top
            invAng=self.shift+self.span
            #                     /
            #                   /
            # involute-angle- /) /
            #                |
            #                |  / - self.top+self.tip
            #       self.top-|
            #                | /
            #              b=|_ -angle to be calculated
            #                |/
            b= invAng-asin(sin(invAng)*self.top/rt)
            dot(rt*cos(phi_start),rt*sin(phi_start))
            phi_start+=b
        x,y=self.draw_arc(phi_start, steps, False)
        #preparation for next cog
        r=self.base*sqrt(1+self.shift**2)
        rot=phi_start+self.span-atan(self.shift+self.span)+atan(self.shift)+2*self.offset

        if self.clearance!=0 and self.bottom == self.pitch-self.D:
            rc=self.bottom-self.clearance
            frac=rc/self.bottom
            #    contact * of involute w. c extension
            #            |\
            #            |_\- angle known from involute(self.shift)
            #            |  \
            #            |   \
            #self.bottom-|    /
            #            |   /- self.bottom-self.clearence
            #            |__/
            #          a=| /- angle to be calculated
            #            |/
            a= -self.shift+asin(sin(self.shift)*self.bottom/rc)
            cx=frac*(x*cos(a)-y*sin(a))
            cy=frac*(x*sin(a)+y*cos(a))
            line(cx,cy, rc*cos(rot-a),rc*sin(rot-a))

        x,y=self.draw_arc(rot, steps)
        if self.tip!=0 and self.bottom == self.pitch-self.D:
            x=frat*(x*cos(b)-y*sin(b))
            y=frat*(x*sin(b)+y*cos(b))
            dot(x,y)
        return x,y

    def draw_arc(self, phi_start, steps, ornt=True):
        sign=(1 if ornt else -1)
        shift=self.shift if ornt else -self.shift-self.span
        phixed=phi_start-(shift-atan(shift))
        #phixed=phi_start+(self.shift+self.span-atan(self.shift+self.span))
        d=self.span/steps
        for i in range(0, steps+1):
            phi= d*i + shift
            sx = cos(phixed+phi) + sin(phixed+phi)*phi
            sy = sin(phixed+phi) - cos(phixed+phi)*phi
            px=sx*self.base
            py=sy*self.base
            dot(px,py)
        return px,py

class Cycloid(Cog):

    def __init__(self, args):
        super(Cycloid, self).__init__(args)
        self.part=(self.m/2+self.spacing)/self.pitch
        top=(self.m/2+self.spacing)/(2*pi)
        low=(self.m/2-self.spacing)/(2*pi)

    def draw(self, phi_start, steps):
        phi_part = self.part+phi_start
        self.draw_arc(phi_start, steps)
        return self.draw_arc(phi_part, steps, False)

    def draw_arc(self, phi_start, steps, ornt=True):
        sign=1 if ornt else -1
        r=(self.m/2+sign*self.spacing)/(2*pi)
        prevx, prevy = cos(phi_start)*self.pitch, sin(phi_start)*self.pitch
        k=1/float(self.N*2)
        base=self.pitch+sign*r
        phi_step=pi/args.N/steps
        for i in range(1,steps+1):
            phi = phi_step*i
            x=base*cos(phi+phi_start) - r*cos(sign*phi_start+base/r*phi)*sign
            y=base*sin(phi+phi_start) - r*sin(sign*phi_start+base/r*phi)
            dot(x,y)
        return x,y

cog=None
if args.involute:
    cog=Involute(args)
else:
    cog=Cycloid(args)

cog.draw_wheel()

if args.outfile:
    if args.outfile != "-":
        f = open(args.outfile, "w")
    else:
        f = sys.stdout

    if not args.e_mode:
        f.write("  0\nSECTION\n")
        f.write("  2\nENTITIES\n")

    for i in range(len(linedata)//2):
        x1 = linedata[i*2 - 2]
        y1 = linedata[i*2 - 1]
        x2 = linedata[i*2 + 0]
        y2 = linedata[i*2 + 1]
        f.write("  0\nLINE\n")
        if args.layername:
            f.write("  8\n{0}\n".format(args.layername))
        f.write(" 10\n{0}\n".format(x1 + args.off_x))
        f.write(" 20\n{0}\n".format(y1 + args.off_y))
        f.write(" 11\n{0}\n".format(x2 + args.off_x))
        f.write(" 21\n{0}\n".format(y2 + args.off_y))

    if not args.e_mode:
        f.write("  0\nENDSEC\n")
        f.write("  0\nEOF\n")

    if args.outfile != "-":
        f.close()

else:
    #import huge dependency only when absolutely necessary!
    import Tkinter
    tk = Tkinter.Tk()
    tk.wm_geometry("400x400")

    canvas = Tkinter.Canvas(tk, width=400, height=400)
    canvas.pack()

    max_xy = 0

    max_xy = max(*linedata)

    scale = 180 / max_xy

    for i in range(len(linedata)//2):
        x1 = linedata[i*2 - 2]
        y1 = linedata[i*2 - 1]
        x2 = linedata[i*2 + 0]
        y2 = linedata[i*2 + 1]
        canvas.create_line(200 + x1*scale, 200 - y1*scale, 200 + x2*scale, 200 - y2*scale)

    x = scale * args.N*args.m/2/pi
    canvas.create_oval(200 - x, 200 - x, 200 + x, 200 + x, outline="green")

    tk.mainloop()

exit(0)
