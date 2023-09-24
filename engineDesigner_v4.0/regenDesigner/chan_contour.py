import numpy as np
import matplotlib.pyplot as plt


def get_chan_contour(jacket):

    rad_inner = jacket.engine.engineContour
    a = np.zeros([len(rad_inner[:, 0]), 1])
    rad_inner = np.concatenate((rad_inner, a), 1)

    b = jacket.engine.engineProps[:-(jacket.start_ind + 1):, 1]

    wall_t = np.concatenate((jacket.wall_t_arr[:, 0], jacket.wall_t_arr[-1, :]), 0)
    channel_h = np.concatenate((jacket.channel_h_arr[:, 0], jacket.channel_h_arr[-1, :]), 0)
    channel_w = np.concatenate((jacket.channel_w_arr[:, 0], jacket.channel_w_arr[-1, :]), 0)
    rad_outer = np.stack((rad_inner[:,0] + wall_t + channel_h, rad_inner[:, 1], rad_inner[:,2]),  axis=1)

    liner_contour = np.concatenate((rad_inner, np.flip(rad_outer, 0), np.array([rad_inner[0, :]])), 0)

    np.savetxt('liner_contour_inner_CR5_PE8_MR2.txt', rad_inner, fmt='%1.6f', delimiter="\t")
    np.savetxt('liner_contour_outer_CR5_PE8_MR2.txt', rad_outer, fmt='%1.6f', delimiter="\t")


    chan_botleft = np.stack((rad_inner[:, 0] + wall_t, rad_inner[:, 1], rad_inner[:, 2] + channel_w / 2), axis=1)
    chan_botright = np.stack((rad_inner[:, 0] + wall_t, rad_inner[:, 1], rad_inner[:, 2] - channel_w / 2), axis=1)

    chan_topleft = np.stack((rad_inner[:, 0] + wall_t + channel_h, rad_inner[:, 1], rad_inner[:, 2] + channel_w / 2), axis=1)
    chan_topright = np.stack((rad_inner[:, 0] + wall_t + channel_h, rad_inner[:, 1], rad_inner[:, 2] - channel_w / 2), axis=1)

    np.savetxt('chan_botleft_0012_60.txt', chan_botleft, fmt='%1.6f', delimiter="\t")
    np.savetxt('chan_botright_0012_60.txt', chan_botright, fmt='%1.6f', delimiter="\t")
    np.savetxt('chan_topleft_0012_60.txt', chan_topleft, fmt='%1.6f', delimiter="\t")
    np.savetxt('chan_topright_0012_60.txt', chan_topright, fmt='%1.6f', delimiter="\t")

    plt.figure()
    plt.plot(liner_contour[:, 0], liner_contour[:, 1])
    #plt.show()


    return