from scipy import optimize

def findPrinciples(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_zx):

    def f(sigma):
        return sigma**3 - (sigma_x + sigma_y + sigma_z) * sigma**2 + (sigma_x * sigma_y + sigma_x * sigma_z + sigma_y * sigma_z - tau_xy**2 - tau_yz**2 - tau_zx**2) * sigma - (sigma_x * sigma_y * sigma_z + 2 * tau_xy * tau_yz * tau_zx - sigma_x * tau_yz**2 - sigma_y * tau_zx**2 - sigma_z * tau_xy**2)

    sigmas = []

    # Find the roots of this function
    for i in range(-1000, 1000):
        # Check every point in the range
        try:
            root = optimize.brentq(f, i, i+1)
            sigmas.append(root) # Add root to list if it exists
        except:
            pass
    sigmas.sort(reverse=True)
    if len(sigmas) == 2:
        sigmas.append(0)
        sigmas.sort(reverse=True)
    elif len(sigmas) == 1:
        sigmas.append([0, 0])
        sigmas.sort(reverse=True)

    return sigmas[0], sigmas[1], sigmas[2]
