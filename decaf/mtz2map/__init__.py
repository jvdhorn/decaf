import sys

from . import mtz2map


def main():
  mtz2map.run(sys.argv[1:])

if __name__ == '__main__':
  main()
