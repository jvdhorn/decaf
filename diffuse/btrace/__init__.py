import sys

from . import btrace


def main():
  btrace.run(sys.argv[1:])

if __name__ == '__main__':
  main()
