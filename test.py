from algebra_stuff import *
import time


def timer(f, *args):
    t0 = time.time()
    res = f(*args)
    t1 = time.time()
    print(f"Executed in {t1-t0} seconds")
    return res

set_global_scope(globals())
R = PolyRing(n=3, make_symbols_global_vars=False)
x, y, z = R.symbols

#J = R.ideal(x**2, x*y**2, x*y*z, x*z**2, y**2*z**2, y*z**3, z**4, y**3-x*z)

f=GroebnerPolynomial.make(x**3+5*x*y**2+2*y*z+z**4+z-9, symbols=R.symbols)

I = R.ideal(x**2, x*y**2, x*y*z, x*z**2, y**2*z**2, y*z**3, z**4, y**3-x*z)
# I = R.ideal(x**2, y**2, z**2)
# I = R.ideal(x, y, z)**2
I1 = ideal(x**2, y, z)
I2 = ideal(x**2, y**2, z**2)

S = R//I
O = R/I
J = I/I**2


s = "x2 x3 x3y x3z x2y x2z xy2 x2y2 x2y2z xy3 xy4 xy3z xy2z xyz x2yz x2yz2 xyz2 xz2 x2z2 x2z3 xz3 y2z2 xy2z2 xy2z3 y3z2 y4z2 y3z3 y2z3 y2z4 y2z5 yz3 xyz3 xyz4 yz4 yz5 yz6 z4 xz4 xz5 z5 z6 z7 y3-xz y4-xyz y5-xy2z y5z-xy2z2 y4z-xyz2 y3z-xz2"


def make_poly(s: str):
    symbols = {'x': x, 'y': y, 'z': z}
    operators = {'+': 1, '-': -1}
    expression = []
    symb = 1
    exp = 1
    p = GroebnerPolynomial.make(1, order=R.order, symbols=R.symbols)
    for c in s:
        if c.isnumeric():
            if symb == 1:
                p *= float(c)
            else:
                exp = int(c)
        elif c in symbols:
            p *= symb**exp
            symb = symbols[c]
            exp = 1
        elif c in operators:
            p *= symb**exp
            expression.append(p)
            p = GroebnerPolynomial.make(operators[c], order=R.order, symbols=R.symbols)
            symb = 1
            exp = 1
    p *= symb**exp
    return sum(expression, start=p)


s = list(map(make_poly, s.split()))
s = R.sort_list(s)


# for Macaulay2: R=QQ[x,y,z]; I=ideal(x^2, x*y^2, x*y*z, x*z^2, y^2*z^2, y*z^3, z^4, y^3-x*z); J=I/I^2; O=R/I;


import subprocess
import threading
import _thread


def execute_commands(commands: List[str]):
    cmd = " && ".join(commands)
    res = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    return res


def interact_with_python():
    try:
        # Start Python interpreter as a subprocess with the -i flag
        python_process = subprocess.Popen(['python', '-i'],
                                         stdin=subprocess.PIPE,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         text=True)

        # Example Python expressions to execute
        expressions = [
            "2 + 2",
            "print('Hello, world!')",
            "import math; math.sqrt(16)"
        ]

        # Execute expressions and capture output
        while True:
            expr = input()
            # Write expression to stdin
            python_process.stdin.write(expr + '\n')
            python_process.stdin.flush()

            # Read output from stdout
            print(repr(python_process.stdout.newlines))
            output = python_process.stdout.readline()
            print("Output:", output.strip())

        # Close stdin to indicate no more input
        python_process.stdin.close()

        # Wait for the process to finish and get its return code
        python_process.wait()
    except Exception as e:
        # Handle errors
        print(f"Error: {e}")


class Rate:
    def __init__(self, hz: int):
        self._hz = hz
        self._period = 1/hz
        self._t = time.time()
    
    def sleep(self):
        t = time.time()-self._t
        time.sleep(max(0, self._period - t))
        self._t = t
    
    def reset(self):
        self._t = time.time()
        

class Macaulay2Prompt:
    def __init__(self):
        self.process = subprocess.Popen(['M2'],
                                        stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        text=True)
        # self._must_read = False
        # self._stdout_line = ""
        self.readlines(0.5)
    
    def close(self):
        self.process.stdin.close()
        self.process.wait()
    
    # def read_thread(self, hz=20):
    #     rate = Rate(hz)
    #     while True:
    #         if self._must_read:
    #             self._stdout_line = self.process.stdout.readline()
    #             self._must_read = False
    #         rate.sleep()
        
    def readline(self, timeout=0.5):
        def _read():
            info[1] = self.process.stdout.readline().strip()
            info[0] = True
        info = [False, None]    # done, line
        threading.Thread(target=_read, daemon=True).start()
        rate = Rate(hz=10)
        t0 = time.time()
        while True:
            if info[0]:
                return info[1]
            if time.time() - t0 > timeout:  # TODO: kill the thread in that case
                return None
            rate.sleep()
    
    def readline2(self, timeout=0.5):
        def interrupt():
            raise TimeoutError
        try:
            timeout_timer = threading.Timer(timeout, _thread.interrupt_main)
            timeout_timer.start()
            line = self.process.stdout.readline().strip()
            return line
        except KeyboardInterrupt:
            return None
        finally:
            timeout_timer.cancel()
    
    def readlines(self, timeout=0.5):
        lines = []
        while (line:=self.readline(timeout)) is not None:
            lines.append(line)
        return lines

    def write(self, expr: str, timeout=0.5):
        self.process.stdin.write(expr + '\n')
        self.process.stdin.flush()
        #output = self.process.stdout.readline()
        #stdout, stderr = self.process.communicate(input=expr + ';\n')
        output = "\n".join(self.readlines(timeout))
        return output
    
    def interact(self):
        while True:
            expr = input("in: ")
            output = self.write(expr)
            print("out:", output.strip())

        
prompt = Macaulay2Prompt()


def test_hom_rank(prompt: Macaulay2Prompt, I: PolyRingIdeal):
    R = I.base
    symbols = R.symbols
    base = "QQ"
    m2_eqs = [f.macaulay2_repr for f in I.gens]
    ring_def = f"R = {base}[{','.join(map(repr, symbols))}]"
    ideal_def = f"I = ideal({','.join(m2_eqs)})"
    degree_comp = "degree Hom(I/I^2, R/I)"
    cmds = [ring_def, ideal_def, degree_comp]
    out = prompt.write(";".join(cmds))
    m2_deg = int(out.strip().split("=")[-1])
    my_deg = hom_rank(I/I**2, R/I)
    print("Macaulay2:", m2_deg)
    print("Computed:", my_deg)
    return m2_deg == my_deg


def random_poly(R: PolyRing, degree: int, with_constant_term=False):
    pass
def random_ideal(R: PolyRing):
    pass


H = DoubleNestedHilbertScheme([2,1,1], R)
