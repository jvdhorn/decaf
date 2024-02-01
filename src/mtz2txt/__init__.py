import sys

from . import mtz2txt


def main():
  mtz2txt.run(sys.argv[1:])

if __name__ == '__main__':
  main()
