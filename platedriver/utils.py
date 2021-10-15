#!/usr/bin/env python3
#platedriver/utils.py


def normalize8(self, I):
        mn = I.min()
        mx = I.max()
        I = ((I - mn)/(mx-mn)) * 255
        return I.astype(np.uint8)

