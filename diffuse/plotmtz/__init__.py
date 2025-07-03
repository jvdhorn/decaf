import sys

from . import plotmtz


def main():
  plotmtz.run(sys.argv[1:])

if __name__ == '__main__':
  main()
