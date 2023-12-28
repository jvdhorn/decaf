import sys

from . import pdbrad


def main():
  pdbrad.run(sys.argv[1:])

if __name__ == '__main__':
  main()
