import sys

from . import pdbdist


def main():
  pdbdist.run(sys.argv[1:])

if __name__ == '__main__':
  main()
