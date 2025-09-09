import numpy as np
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Take a input numpy array file and output its shape"
    )
    p.add_argument(
        'input_numpy_file',
        help=".npy file need to be examinated",
    )
    return p.parse_args()

def main():
    args = parse_args()
    arr = np.load(args.input_numpy_file)
    print(arr.shape)


if __name__ == '__main__':
    main()