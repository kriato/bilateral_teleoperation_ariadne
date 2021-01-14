mod = 9
a = 2
c = 0
seed = 1

# https://en.m.wikipedia.org/wiki/Linear_congruential_generator
def lcg():
    global mod, a, c, seed

    # seed is the random variable
    seed = (a * seed + c) % mod
    mod += 1
    a += 1
    c += 1

def lcg_nomod():
    global mod, a, c, seed

    # seed is the random variable
    # modulus is implemented as a - b * int(a/b)
    var = a * seed + c
    seed = var - mod * int(var / mod)
    mod += 1
    a += 1
    c += 1

it = 10
success = []
loss = []
all_it = []
for _ in range(it):
    # lcg()
    lcg_nomod()
    print(seed)
    all_it.append(seed)

    if seed % 3 == 0:
        success.append(1)
        loss.append(0)
    else:
        loss.append(1)
        success.append(0)
   
import matplotlib.pyplot as plt
import numpy as np

plt.plot(all_it)
plt.show()

print(f'Success: {np.count_nonzero(success)/ float(it) * 100}%')
print(f'Loss: {np.count_nonzero(loss)/ float(it) * 100}%')