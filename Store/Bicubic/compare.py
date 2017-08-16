import numpy as np

_logTe = np.log10(50)
_logNe = np.log10(8e18)
_k = 0

u = np.array(
[[-18.601440, -18.751220, -18.968840, -19.133310],
[-18.565540, -18.715990, -18.934000, -19.099460],
[-18.530450, -18.680260, -18.899830, -19.066350],
[-18.485350, -18.634280, -18.855230, -19.022510]])

print("log_coeff subgrid")
print(u)

# Switch to [x, y] indexing
u = np.transpose(u)
u = np.flip(u,1)

logTe = np.array([1.301230, 1.477320, 1.699170, 1.845300]);
logNe = np.array([18.301030, 18.698970, 19.000000, 19.301030]);

print("Te: ",logTe)
print("Ne: ",logNe)

f_sub = np.zeros((4,4));
f__ = np.zeros((2,2));
f_t = np.zeros((2,2));
f_n = np.zeros((2,2));
ftn = np.zeros((2,2));

# Have to swap order, which is a bit confusing
shift = 1;

for x in range(1,3):
	for y in range(1,3):
		i = y-shift
		j = x-shift

		f__[i,j] = u[x, y]

		f_t[i,j] = (u[(x+1), y]     - u[(x-1), y])     / (logTe[(x+1)]-logTe[(x-1)])
		
		# print("f_t({}, {}) {:+2.5f} = ({:+2.5f} - {:+2.5f}) / ({:+2.5f}-{:+2.5f})\n".format(i,j,
		# f_t[i,j],
		# u[(x+1), y]     , u[(x-1), y]     , logTe[(x+1)],logTe[(x-1)]))

		f_n[i,j] = (u[x,     (y+1)] - u[x,     (y-1)]) / (logNe[(y+1)]-logNe[(y-1)])
		
		# print("f_n({}, {}) {:+2.5f} = ({:+2.5f} - {:+2.5f}) / ({:+2.5f}-{:+2.5f})\n".format(i,j,
		# f_n[i,j],
		# u[x,     (y+1)], u[x,     (y-1)], logNe[(y+1)],logNe[(y-1)]))

		ftn[i,j] = (u[x+1, y+1] - u[x+1, y-1] - u[x-1, y+1] + u[x-1, y-1]) / ((logTe[(x+1)]-logTe[(x-1)])*(logNe[(y+1)]-logNe[(y-1)]))

f_sub[0,0] = f__[0,0]
f_sub[0,1] = f__[1,0]
f_sub[1,0] = f__[0,1]
f_sub[1,1] = f__[1,1]

f_sub[0+2,0] = f_t[0,0]
f_sub[0+2,1] = f_t[1,0]
f_sub[1+2,0] = f_t[0,1]
f_sub[1+2,1] = f_t[1,1]

f_sub[0,0+2] = f_n[0,0]
f_sub[0,1+2] = f_n[1,0]
f_sub[1,0+2] = f_n[0,1]
f_sub[1,1+2] = f_n[1,1]

f_sub[0+2,0+2] = ftn[0,0]
f_sub[0+2,1+2] = ftn[1,0]
f_sub[1+2,0+2] = ftn[0,1]
f_sub[1+2,1+2] = ftn[1,1]

print("f_sub")
print(f_sub)

prematrix = np.array(
[[+1, +0, +0, +0],
[+0, +0, +1, +0],
[-3, +3, -2, -1],
[+2, -2, +1, +1]]);

postmatrix = np.array(
[[+1, +0, -3, +2],
[+0, +0, +3, -2],
[+0, +1, -2, +1],
[+0, +0, -1, +1]]);

alpha_sub = np.dot(np.dot(prematrix,f_sub),postmatrix)

print("alpha_sub")
print(alpha_sub)

lowTe = np.searchsorted(logTe,_logTe,side='left')-1
lowNe = np.searchsorted(logNe,_logNe,side='left')-1

print("lowTe: ", lowTe, " -> ({:.4e} < {:.4e} < {:.4e})".format(logTe[lowTe],_logTe,logTe[lowTe+1]))
print("lowNe: ", lowNe, " -> ({:.4e} < {:.4e} < {:.4e})".format(logNe[lowNe],_logNe,logNe[lowNe+1]))

x = (_logTe - logTe[lowTe])/(logTe[lowTe+1] - logTe[lowTe])
y = (_logNe - logNe[lowNe])/(logNe[lowNe+1] - logNe[lowNe])

x_vector = np.array([1, x, x**2, x**3])
y_vector = np.array([1, y, y**2, y**3])
# y_vector = y_vector.reshape((1,4))

print("x_vector")
print(x_vector)
print("y_vector")
print(y_vector)

return_value = np.dot(np.dot(x_vector,alpha_sub),y_vector)

print("return_value")
print(return_value)



