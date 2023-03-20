----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date: 11/17/2022 09:37:33 AM
-- Design Name: 
-- Module Name: add_node_v3 - Behavioral
-- Project Name: 
-- Target Devices: 
-- Tool Versions: 
-- Description: 
-- 
-- Dependencies: 
-- 
-- Revision:
-- Revision 0.01 - File Created
-- Additional Comments:
-- 
----------------------------------------------------------------------------------


library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx leaf cells in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity add_node_v3 is
    generic (
        w_x_i : integer;
        w_y_i : integer;
        w_o : integer;
        s_x_i : integer;
        s_y_i : integer;
        s_o : integer;
        copy_sign_x_i : boolean := false;
        copy_sign_y_i : boolean := false;
        sub : boolean
    );
    port (
        x_i : in signed(w_x_i-1 downto 0);
        y_i : in signed(w_y_i-1 downto 0);
        z_o : out signed(w_o-1 downto 0)
    );
end add_node_v3;

architecture add_node_v3 of add_node_v3 is

impure function get_result_width(s_x : integer; s_y : integer) return integer is
    variable w_temp : integer;
begin
    if s_x + w_x_i >= s_y + w_y_i then
        w_temp := s_x + w_x_i;
    else
        w_temp := s_y + w_y_i;
    end if;
    if s_o > 0 then
        w_temp := w_temp + s_o;
    end if;
    if w_temp >= w_o then
        return w_temp;
    else
        return w_o;
    end if;
end function;

constant result_width : integer := get_result_width(s_x_i, s_y_i);

signal x_ext : signed(result_width-1 downto 0);
signal y_ext : signed(result_width-1 downto 0);

signal x_shifted : signed(result_width-1 downto 0);
signal y_shifted : signed(result_width-1 downto 0);

signal z_res : signed(result_width-1 downto 0);
signal z_sign_copy : signed(result_width-1 downto 0);
signal z_shifted : signed(result_width-1 downto 0);

begin
    -- resize inputs
    x_ext <= resize(x_i, result_width);
    y_ext <= resize(y_i, result_width);
    
    -- shift x
    gen_shift_x : if s_x_i > 0 generate
        x_shifted <= shift_left(x_ext, s_x_i);
    end generate;
    gen_no_shift_x : if s_x_i = 0 generate
        x_shifted <= x_ext;
    end generate;
    
    -- shift y
    gen_shift_y : if s_y_i > 0 generate
        y_shifted <= shift_left(y_ext, s_y_i);
    end generate;
    gen_no_shift_y : if s_y_i = 0 generate
        y_shifted <= y_ext;
    end generate;

    -- calc output
    gen_sub : if sub = true generate
        z_res <= x_shifted - y_shifted;
    end generate;
    gen_add : if sub = false generate
        z_res <= x_shifted + y_shifted;
    end generate;
    
    -- handle sign
    gen_sign_copy_x : if copy_sign_x_i generate
        z_o(w_o-1) <= x_i(w_x_i-1);
    end generate;
    gen_sign_copy_y : if copy_sign_y_i generate
        z_o(w_o-1) <= y_i(w_y_i-1);
    end generate;
    gen_no_sign_copy : if not (copy_sign_x_i or copy_sign_y_i) generate
        z_o(w_o-1) <= z_res(w_o-1);
    end generate;
    
    -- handle remaining bits depending on post-add-shift
    z_o(w_o-2 downto 0) <= z_res(w_o-2+s_o downto s_o);

end add_node_v3;
