{ pkgs, lib, config, inputs, ... }:

{
  languages.python = {
    enable = true;
    poetry.enable = true;
    poetry.install.enable = true;
  };
}
