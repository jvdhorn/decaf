import sys

from . import subadp


def main():
  subadp.run(sys.argv[1:])

if __name__ == '__main__':
  main()
