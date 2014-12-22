/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package main;

import canvas.DisplayFrame;
import canvas.InputFrame;
import control.Plotter;

/**
 *
 * @author edsonlan
 */
public class TestNativeCanvas {

    public static void main(String args[]) {

        System.loadLibrary("rpn");

        Numerics.initNumerics();

        DisplayFrame frame = new DisplayFrame("Frame 1",1);
        
//        DisplayFrame frame2 = new DisplayFrame("Frame 2", 2);        
//        
//        DisplayFrame frame3 = new DisplayFrame("Frame 3", 3);        

        
        
        
        
        Plotter plotter = new Plotter();
        plotter.add(frame);
//        plotter.add(frame2);
//        plotter.add(frame3);

        InputFrame input = new InputFrame("Input",plotter);
//
        input.setVisible(true);

        frame.setVisible(true);
        
//        frame2.setVisible(true);
//        
//        frame3.setVisible(true);
        
        

    }

}
