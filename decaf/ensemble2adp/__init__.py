import sys

from . import ensemble2adp


def main():
  ensemble2adp.run(sys.argv[1:])

if __name__ == '__main__':
  main()
