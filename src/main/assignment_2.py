import numpy as np

# neville's method for lagrangian interpolation
def neville_method(points, x):
  numPoints = len(points[0])
  table = np.zeros((numPoints,numPoints))
  table[0] = points[1]

  for i in range(1,numPoints):
    for j in range(i,numPoints):
      table[i][j] = ((x-points[0][j-i])*table[i-1][j]-(x-points[0][j])*table[i-1][j-1]) /(points[0][j]-points[0][j-i])
  print(table[numPoints-1][numPoints-1], end="\n\n")


# newton's forward method for lagrangian interpolation
def newton_method(points,x):
  numPoints = len(points[0])
  table = np.zeros((numPoints,numPoints))
  table[0] = points[1]

  for i in range(1,numPoints):
    for j in range(i,numPoints):
      table[i][j] = (table[i-1][j] - table[i-1][j-1]) / (points[0][j]-points[0][j-i])
  
  print("[", end="")
  for i in range(1,numPoints):
    print(table[i][i], end = '')
    if i < numPoints-1:
      print(", ", end = '')
  print("]", end ='\n\n')

  y = 0
  for i in range(0, len(points[0])):
    product = table[i][i]
    for j in range(0,i):
      product *= (x - points[0][j])
    y += product

  print(y, end='\n\n')

# hermite approximation with divided differences table
def hermite_method(points):
  numPoints = len(points[0])
  table = np.zeros((numPoints*2,numPoints*2))

  for i in range(0,2):
    for j in range(0, numPoints):
      table[2*j][i] = points[i][j]
      table[2*j+1][i] = points[i][j]
  
  for j in range(0,numPoints):
    table[2*j+1][2] = points[2][j]
  for j in range(2,numPoints*2,2):
    table[j][2] = (table[j][1]-table[j-1][1])/(table[j][0]-table[j-1][0])

  for i in range(3, 2*numPoints):
    for j in range(i-1, 2*numPoints):
      table[j][i] = (table[j][i-1]-table[j-1][i-1])/(table[j][0]-table[j+1-i][0])

  print(table, end='\n\n')

# cubic spline interpolation
def cubic_spline(points):
  n = len(points[0])
  h = np.empty([n-1])
  for i in range(0,n-1):
    h[i] = points[0,i+1] - points[0,i]
  
  A = np.zeros((n,n))
  A[0,0] = 1
  A[n-1,n-1] = 1

  for i in range(1,n-1):
    A[i,i-1] = h[i-1]
    A[i,i] = 2*(h[i-1]+h[i])
    A[i,i+1] = h[i]
  
  print(A, end='\n\n')

  b = np.zeros([n])
  for i in range(1,n-1):
    b[i] = 3/h[i]*(points[1,i+1]-points[1,i])-3/h[i-1]*(points[1,i]-points[1,i-1])

  print(b, end = '\n\n')

  c = np.matmul(np.linalg.inv(A),b)

  print(c, end ='\n\n')

np.set_printoptions(precision=7, suppress=True, linewidth=100)

neville_method([[3.6,3.8,3.9],[1.675,1.436,1.318]],3.7)

newton_method([[7.2,7.4,7.5,7.6],[23.5492,25.3913,26.8224,27.4589]],7.3)

hermite_method([[3.6,3.8,3.9],[1.675,1.436,1.318],[-1.195,-1.188,-1.182]])

cubic_spline(np.array([[2,5,8,10],[3,5,7,9]]))