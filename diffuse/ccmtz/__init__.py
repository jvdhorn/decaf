import sys

from . import ccmtz


def main():
  ccmtz.run(sys.argv[1:])

if __name__ == '__main__':
  main()
