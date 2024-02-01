import sys

from . import rmbragg


def main():
  rmbragg.run(sys.argv[1:])

if __name__ == '__main__':
  main()
