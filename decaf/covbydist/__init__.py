import sys

from . import covbydist


def main():
  covbydist.run(sys.argv[1:])

if __name__ == '__main__':
  main()
