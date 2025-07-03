import sys

from . import patterson


def main():
  patterson.run(sys.argv[1:])

if __name__ == '__main__':
  main()
