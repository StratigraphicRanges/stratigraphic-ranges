package sr.treeannotator;

import java.util.Optional;

import beast.base.core.Log;
import beastfx.app.beauti.ThemeProvider;
import beastfx.app.inputeditor.BEASTObjectDialog;
import beastfx.app.inputeditor.BEASTObjectPanel;
import beastfx.app.inputeditor.BeautiConfig;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.Console;
import beastfx.app.util.Utils;
import javafx.scene.control.ButtonType;

/**
 * GUI front end for {@link SRTreeAnnotator}, modelled on BEAST's TAGui.
 *
 * Shows the same input dialog as beastfx.app.tools.Application would, but
 * runs the annotation on a worker thread instead of the JavaFX application
 * thread. The console window therefore repaints while the analysis runs and
 * progress output appears live instead of all at once at the end.
 */
public class SRTreeAnnotatorGui extends Console {

    private static final String TITLE = "SR Tree Annotator";

    public static void launchGui(String[] args) {
        launch(SRTreeAnnotatorGui.class, args);
    }

    @Override
    protected void createDialog() {
        if (Utils.isMac()) {
            Utils.loadUIManager();
        }

        BeautiDoc doc = new BeautiDoc();
        doc.beautiConfig = new BeautiConfig();
        doc.beautiConfig.initAndValidate();

        SRTreeAnnotator annotator = new SRTreeAnnotator();
        annotator.setID(TITLE);

        BEASTObjectPanel panel = new BEASTObjectPanel(annotator, annotator.getClass(), doc);
        BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
        dialog.setResizable(true);
        dialog.getDialogPane().setPrefSize(1024, 768);
        dialog.getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        ThemeProvider.loadStyleSheet(dialog.getDialogPane().getScene());

        Optional<ButtonType> option = dialog.showAndWait();
        if (!option.isPresent() || !option.get().equals(ButtonType.OK)) {
            Log.warning("Quitting " + TITLE);
            System.exit(0);
            return;
        }
        dialog.accept(annotator, doc);

        Thread worker = new Thread(() -> {
            long start = System.currentTimeMillis();
            try {
                annotator.initAndValidate();
                annotator.run();
            } catch (Exception e) {
                Log.err.println(TITLE + " failed: " + e.getMessage());
                e.printStackTrace();
            }
            long end = System.currentTimeMillis();
            Log.warning(TITLE + " done in " + (end - start) / 1000
                    + " seconds. Close window to finish application.");
        }, "SRTreeAnnotator-worker");
        worker.start();
    }
}
