import sys

from . import map2mtz


def main():
  map2mtz.run(sys.argv[1:])

if __name__ == '__main__':
  main()
