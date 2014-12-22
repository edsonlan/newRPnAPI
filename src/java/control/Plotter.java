/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package control;

import canvas.DisplayFrame;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JFrame;

/**
 *
 * @author edsonlan
 */
public class Plotter {

    private List<JFrame> displayFramesList_;

    public Plotter() {

        displayFramesList_ = new ArrayList<JFrame>();

    }

    public void add(JFrame frame) {
        displayFramesList_.add(frame);
    }

    public void update() {

        for (JFrame jFrame : displayFramesList_) {

            DisplayFrame displayFrame = (DisplayFrame) jFrame;

            displayFrame.getCanvas().paint(displayFrame.getCanvas().getGraphics());

        }

    }
    
    
    
    
    public void cleanFrames(){
        
        for (JFrame jFrame : displayFramesList_) {
            
            jFrame.dispose();
            
        }
        
        
        
    }

    public native void cleanList();
    
    public native void plot(String curveName, int xDC,int yDC);

}
