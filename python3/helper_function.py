# Helper functions.
import matplotlib.pyplot as plt
import numpy as np
import shutil
from scipy import stats


###############
# MULTITHREAD #
###############

def worker(obj_arg):
    obj = obj_arg[0]
    arg = obj_arg[1:]
    print(arg)
    return obj.simul(*arg)



##############
# STATISTICS #
##############

def t_test_unpaired(a,b, tail = 1):
    """
    Unpaired t-test (one or two tailed),
    If one tail, Check if a > to b.
    the 2 tails version is equivalent to stats.ttest_ind(a,b).
    """
    # Check if Eest[a] > Eest[b]
    N = len(a)
    # For unbiased max likelihood estimate we have to divide the var by N-1,
    # and therefore the parameter ddof = 1
    var_a = a.var(ddof=1)
    var_b = b.var(ddof=1)
    print(var_a, var_b)
    #std deviation
    s = np.sqrt((var_a + var_b)/2)
    print("standard deviation " + str(s))
    ## Calculate the t-statistics
    t = (a.mean() - b.mean())/(s*np.sqrt(2/float(N)))
    ## Compare with the critical t-value
    #Degrees of freedom
    df = 2*N - 2
    #p-value after comparison with the t
    p = 1 - stats.t.cdf(t,df=df)
    print(("t = " + str(t)))
    if tail == 1:
        print(("p = " + str(p)))
    if tail == 2:
        print(("p = " + str(2*p)))

    return t,p

def t_test_unpaired(a,b, tail = 1):
    """
    Update  1: This version accept different length for the 2 samples.
    Update 2: Accept list or vector as input.
    Unpaired t-test (one or two tailed),
    If one tail, Check if a > to b.
    the 2 tails version is equivalent to stats.ttest_ind(a,b).
    """
    # Check if Eest[a] > Eest[b]
    Na = len(a)
    Nb = len(b)
    # For unbiased max likelihood estimate we have to divide the var by N-1,
    # and therefore the parameter ddof = 1
    try:
        var_a = a.var(ddof=1)
        var_b = b.var(ddof=1)
    except:
        a,b = np.array(a), np.array(b)
        var_a = a.var(ddof=1)
        var_b = b.var(ddof=1)


    print(var_a, var_b)
    #std deviation
    s = np.sqrt(((Na-1)*var_a + (Nb-1)*var_b)/(Na + Nb -2))
    print("standard deviation " + str(s))
    ## Standard error of the difference of the mean.
    se_dm = (s*np.sqrt(1/float(Na) + 1/float(Nb)))
    ## Calculate the t-statistics
    t = (a.mean() - b.mean())/se_dm
    ## Compare with the critical t-value
    #Degrees of freedom
    df = Na + Nb - 2
    #p-value after comparison with the t
    p = 1 - stats.t.cdf(t,df=df)
    print(("t = " + str(t)))
    if tail == 1:
        print(("p = " + str(p)))
    if tail == 2:
        print(("p = " + str(2*p)))

    return t,p



######################
# NUMPY COMPUTATIONS #
######################

def permut(N, r):
    """
    Return an array of r*N random numbers between 0 and N (no duplicates).
    """
    N_in_permut = int(N * r)
    N_permut = np.random.permutation(list(range(N)))
    return N_permut[:N_in_permut]


def rolling_average(X, w, dim = 0, method = 'pd'):
    if method == 'pd':
        # Here, we are using Pandas method because fast and easy to set up.
        # However it is depreciated and will be deleted from future release
        # of pandas...
        X_av = pd.rolling_mean(X, w)
        X_av = X_av[w-1:]
    if method == 'naive':
        # First method tried: not efficient at all.
        lX = X.shape[dim]
        tmp = [X[i:i+w].mean(dim) for i in range(lX - w + 1)]
        X_av = np.array(tmp)
    return X_av

# def fast_rolling_average(X, w, dim = 0):
#     pd.rolling_mean(X, N)[N-1:]

def equal_vec(v1, v2):
    """
    Check if all elements of 2 vectors (of the same length) are equal.
    """
    return (v1 == v2).sum() == v1.size

def distance(position1, position2):
    """
    Return the distance between two positions.
    Can be used with 2D positions, e.g. (x0, y0) and (x1, y1)
    or 3D positions, e.g. (x0, y0, z0) and (x1, y1, z1)
    """
    tmp = 0
    for i in range(len(position1)):
        tmp += (position1[i] - position2[i])**2
    return np.sqrt(tmp)

# Mean Square Error.
def MSE(yest, y):
    # Possible to divide that by 2 (for convenience when derivate).
    tmp = (yest - y)**2
    return tmp.mean()


def interval(boolvec):
    """
    Gives the interval between each occurence of '1' or 'True' in a boolvec.
    """
    ### index of the '1'.
    idx1 = np.where(boolvec)[0]
    ### compute intervals.
    interval = [idx1[i+1] - idx1[i] for i in range(len(idx1) - 1)]
    return interval

def in_tupleRange(x, tupleRange):
    """
    Return True if tupleRange[0] <= x <= tupleRange[1]
    """
    return (tupleRange[0] <= x) and (x <= tupleRange[1])

def in_start_dt(x, start, dt):
    """
    Return True if start <= x <= start + dt
    """
    return (start <= x) and (x <= start + dt)


###########
# GENERAL #
###########


def stringDate():
    """
    Return a string = "yyyymdhms"
    """
    import time
    strDate = ""
    for c in time.localtime()[:-3]:
        strDate += str(c)
    return strDate

def mean_list(l):
    return sum(l) / float(len(l))

def add_dim(v):
    """
    Transform a shape (l,) into (l,1)
    ! Careful, with this code also change (,l) into (l,1) !
    = Identity if shape[1] => 1
    """
    try:
        # If there is no error it means we don't need to change the shape.
        d = v.shape[1]
        return v
    except:
        return v.reshape((len(v), 1))

def number_from_list(l):
    """
    Return a number from a list (e.g. [1,2,3] -> 123).
    """
    ll = len(l)
    num = 0
    for i in range(ll):
        num += l[i] * pow(10, ll - i -1)
    return num

def intList_from_string(string):
    intList = []
    for e in string:
        try:
            intList.append(int(e))
        except:
            pass
    return intList

# def fill_matrix_from_lists_if_in_dict(row_list, column_list, dic_of_dic):
#     """
#     For each element r_i of row_list, if dic_of_dic[r_i] exists, it look for
#     each element r_j of column_list if dic_of_dic[r_i][c_j] exists. If it does
#     we stock it in output_matrix[i][j].
#     """
#     Nrow = len(row_list)
#     Ncolumn = len(column_list)
#     out_mat = np.zeros((Nrow, Ncolumn))
#     for i in range(Nrow):
#         r_i = row_list[i]
#         print r_i
#         try:
#             dic_i = dic_of_dic[r_i]
#             # If dic_i exists, let's continue, otherwise let's try dic_i+1
#             print dic_i
#             for j in range(Ncolumn):
#                 r_j = column_list[j]
#                 try:
#                     out_mat[i][j] = dic_i[r_j]
#                     # if dic_i[r_j] exists, let's continue, otherwise let's
#                     # try dic_i[r_j+1]
#                     print dic_i[r_j]
#                 except:
#                     #print "Houra1"
#                     #pass
#                     out_mat[i][j] = np.inf
#         except:
#             # pass
#             #print "Noo"
#             out_mat[i].fill(np.inf)
#     return out_mat



### Visualization.
# def plot_values(values, label = "", show_flag = True):
#     l = len(values)
#     t = range(l)
#     plt.plot(t, values, label = label)
#     if show_flag:
#         plt.show()

def plot_values(values, label = "", linestyle = "-",
                color = None, show_flag = True):
    l = len(values)
    t = list(range(l))
    if color == None:
        plt.plot(t, values, linestyle = linestyle, label = label)
    else:
        plt.plot(t, values, linestyle = linestyle, color = color,
                 label = label)
    if label != "":
        plt.legend()
    if show_flag:
        plt.show()



def plot_scatter(values, label = "", show_flag = True):
    l = len(values)
    t = list(range(l))
    plt.scatter(t, values, label = label)
    if show_flag:
        plt.show()




def scatter_from_mat(X, size_marker = 10.):
    plt.close()
    Ntimesteps = X.shape[0]
    Nneuron = X.shape[1]
    size_simu = float(Ntimesteps)
    for t in range(Ntimesteps):
        who_fire = np.where(X[t])[0]
        t_vec = np.ones(len(who_fire)) * t
        plt.plot(t_vec, who_fire, marker='o', linestyle='none',
                 color='b', markersize = 100. * size_marker / size_simu)
    plt.show()


#####################
# FILE MANIPULATION #
#####################

def create_or_replace_directory(path):
    """
    Create a new directory (and eventually overwrite old directory with same name).
    ! Careful ! this will erase a directory and its contents if there was
    one named 'path'
    """
    if os.path.exists(path):
        shutil.rmtree(path)

    os.makedirs(path)


### Write a dictionary into a file.
def write_dictionary(dic, path):
    with open(path, 'w') as file:
        file.write(dic)
