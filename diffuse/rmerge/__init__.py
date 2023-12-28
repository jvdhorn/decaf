import sys

from . import rmerge


def main():
  rmerge.run(sys.argv[1:])

if __name__ == '__main__':
  main()
