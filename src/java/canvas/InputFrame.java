/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package canvas;

import control.Plotter;
import java.awt.HeadlessException;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import javax.swing.JFrame;
import main.Numerics;

/**
 *
 * @author edsonlan
 */
public class InputFrame extends JFrame implements WindowListener{
    
    
    private InputCanvas canvas_;
    private Plotter plotter_;

    public InputFrame(String title,Plotter plotter) throws HeadlessException {
        super(title);
        addWindowListener(this);
        plotter_=plotter;
        canvas_= new InputCanvas(plotter);
        getContentPane().add(canvas_);
        setSize(400,400);
    }

    @Override
    public void windowOpened(WindowEvent e) {

    }

    @Override
    public void windowClosing(WindowEvent e) {

        dispose();
    }

    @Override
    public void windowClosed(WindowEvent e) {
        Numerics.clean();
        plotter_.cleanFrames();

    }

    @Override
    public void windowIconified(WindowEvent e) {

    }

    @Override
    public void windowDeiconified(WindowEvent e) {

    }

    @Override
    public void windowActivated(WindowEvent e) {

    }

    @Override
    public void windowDeactivated(WindowEvent e) {

    }

    
    
}
