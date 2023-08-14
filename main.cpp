#include "include/window.h"

int main(int argc, char *argv[])
{
  Mesh mesh("../data/spot.obj");
  Window window(mesh);
  window.viewer.launch();

  return 0;
}
